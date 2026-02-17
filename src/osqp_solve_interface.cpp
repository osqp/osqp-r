

#include "osqp.h"
#include "types.h"
#include "algebra_impl.h"
#include <Rcpp.h>
#include <memory>
#include <algorithm>
#include <vector>

using namespace Rcpp;

void extractMatrixData(const Rcpp::S4& mat, std::vector<OSQPInt>& iout, std::vector<OSQPInt>& pout, std::vector<OSQPFloat>& xout);
void translateSettings(OSQPSettings* settings, const Rcpp::List& pars);
void mycleanup(OSQPSolver* x);
S4 toDgCMat(OSQPCscMatrix* inmat);

// predicate for testing if a value is below -OSQP_INFTY
bool below_osqp_neg_inf(OSQPFloat x) {
  return (x < -OSQP_INFTY);
}

// predicate for testing if a value is above OSQP_INFTY
bool above_osqp_inf(OSQPFloat x) {
  return (x > OSQP_INFTY);
}

// [[Rcpp::export]]
SEXP osqpSetup(const S4& P, const NumericVector& q, const S4& A, const NumericVector& l, const NumericVector& u, const List& pars)
{
  IntegerVector dimP = P.slot("Dim");
  IntegerVector dimA = A.slot("Dim");
  OSQPInt n = dimP[0];
  OSQPInt m = dimA[0];
  if (n != dimP[1] || n != dimA[1]) stop("Dimension mismatch in P or A");

  std::vector<OSQPInt> A_i, A_p, P_i, P_p;
  std::vector<OSQPFloat> A_x, P_x, qvec(q.size()), lvec(l.size()), uvec(u.size());

  extractMatrixData(P, P_i, P_p, P_x);
  extractMatrixData(A, A_i, A_p, A_x);

  std::copy(q.begin(), q.end(), qvec.begin());
  std::copy(l.begin(), l.end(), lvec.begin());
  std::copy(u.begin(), u.end(), uvec.begin());

  // Threshold lvec to range [-OSQP_INFTY, OSQP_INFTY]
  std::replace_if(lvec.begin(), lvec.end(), below_osqp_neg_inf, -OSQP_INFTY);
  std::replace_if(lvec.begin(), lvec.end(), above_osqp_inf, OSQP_INFTY);

  // Threshold uvec to range [-OSQP_INFTY, OSQP_INFTY]
  std::replace_if(uvec.begin(), uvec.end(), below_osqp_neg_inf, -OSQP_INFTY);
  std::replace_if(uvec.begin(), uvec.end(), above_osqp_inf, OSQP_INFTY);

  // Create CSC matrices on the stack, referencing our vectors' data
  OSQPCscMatrix Pcsc, Acsc;
  OSQPCscMatrix_set_data(&Pcsc, n, n, P_x.size(), P_x.data(), P_i.data(), P_p.data());
  OSQPCscMatrix_set_data(&Acsc, m, n, A_x.size(), A_x.data(), A_i.data(), A_p.data());

  // Create settings with defaults, then apply user parameters
  OSQPSettings settings;
  osqp_set_default_settings(&settings);
  if (pars.size())
    translateSettings(&settings, pars);

  // Setup solver
  OSQPSolver* solverp = nullptr;
  OSQPInt exitflag = osqp_setup(&solverp, &Pcsc, qvec.data(), &Acsc, lvec.data(), uvec.data(), m, n, &settings);

  if (exitflag != 0 || solverp == nullptr) {
    if (solverp != nullptr) osqp_cleanup(solverp);
    stop("OSQP setup failed with error code %d: %s", exitflag, osqp_error_message(exitflag));
  }

  Rcpp::XPtr<OSQPSolver, Rcpp::PreserveStorage, mycleanup> solver(solverp);
  return solver;
}


// [[Rcpp::export]]
List osqpSolve(SEXP workPtr)
{
  auto solver = as<Rcpp::XPtr<OSQPSolver, Rcpp::PreserveStorage, mycleanup>>(workPtr);

  OSQPInt m_dim, n_dim;
  osqp_get_dimensions(solver, &m_dim, &n_dim);

  // Two-level check: osqp_solve returns 0 on success, nonzero on C-level error
  OSQPInt exitflag = osqp_solve(solver);

  if (exitflag != 0) {
    stop("OSQP solve failed with error code %d: %s", exitflag, osqp_error_message(exitflag));
  }

  std::string status = solver->info->status;
  List info = List::create(
    _("iter") = solver->info->iter,
    _("status") = status,
    _("status_val") = solver->info->status_val,
    _("status_polish") = solver->info->status_polish,
    _("obj_val") = solver->info->obj_val,
    _("dual_obj_val") = solver->info->dual_obj_val,
    _("prim_res") = solver->info->prim_res,
    _("dual_res") = solver->info->dual_res,
    _("duality_gap") = solver->info->duality_gap,
    _("setup_time") = solver->info->setup_time,
    _("solve_time") = solver->info->solve_time,
    _("update_time") = solver->info->update_time,
    _("polish_time") = solver->info->polish_time,
    _("run_time") = solver->info->run_time,
    _("rho_estimate") = solver->info->rho_estimate,
    _("rho_updates") = solver->info->rho_updates);

  NumericVector x(solver->solution->x, solver->solution->x + n_dim);
  NumericVector y(solver->solution->y, solver->solution->y + m_dim);
  NumericVector prim_inf_cert(solver->solution->prim_inf_cert, solver->solution->prim_inf_cert + m_dim);
  NumericVector dual_inf_cert(solver->solution->dual_inf_cert, solver->solution->dual_inf_cert + n_dim);

  List resl = List::create(
    _("x") = x,
    _("y") = y,
    _("prim_inf_cert") = prim_inf_cert,
    _("dual_inf_cert") = dual_inf_cert,
    _("info") = info);

  return resl;
}

// [[Rcpp::export]]
List osqpGetParams(SEXP workPtr)
{
  auto solver = as<Rcpp::XPtr<OSQPSolver, Rcpp::PreserveStorage, mycleanup>>(workPtr);

  IntegerVector linsys;
  if (solver->settings->linsys_solver == OSQP_DIRECT_SOLVER)
    linsys = IntegerVector::create(_("OSQP_DIRECT_SOLVER") = solver->settings->linsys_solver);
  else if (solver->settings->linsys_solver == OSQP_INDIRECT_SOLVER)
    linsys = IntegerVector::create(_("OSQP_INDIRECT_SOLVER") = solver->settings->linsys_solver);
  else
    linsys = IntegerVector::create(_("OSQP_UNKNOWN_SOLVER") = solver->settings->linsys_solver);

  List res = List::create(
    _("rho") = solver->settings->rho,
    _("sigma") = solver->settings->sigma,
    _("max_iter") = solver->settings->max_iter,
    _("eps_abs") = solver->settings->eps_abs,
    _("eps_rel") = solver->settings->eps_rel,
    _("eps_prim_inf") = solver->settings->eps_prim_inf,
    _("eps_dual_inf") = solver->settings->eps_dual_inf,
    _("alpha") = solver->settings->alpha,
    _("linsys_solver") = linsys,
    _("delta") = solver->settings->delta,
    _("polishing") = (bool)solver->settings->polishing,
    _("polish_refine_iter") = solver->settings->polish_refine_iter,
    _("verbose") = (bool)solver->settings->verbose,
    _("scaled_termination") = (bool)solver->settings->scaled_termination,
    _("check_termination") = solver->settings->check_termination,
    _("warm_starting") = (bool)solver->settings->warm_starting,
    _("scaling") = solver->settings->scaling,
    _("adaptive_rho") = solver->settings->adaptive_rho,
    _("adaptive_rho_interval") = solver->settings->adaptive_rho_interval,
    _("adaptive_rho_tolerance") = solver->settings->adaptive_rho_tolerance);

  // 20 is the limit for the List::create construct, so adding others below
  res.push_back(solver->settings->adaptive_rho_fraction, "adaptive_rho_fraction");
  res.push_back(solver->settings->time_limit, "time_limit");
  res.push_back((bool)solver->settings->rho_is_vec, "rho_is_vec");
  res.push_back(solver->settings->cg_max_iter, "cg_max_iter");
  res.push_back(solver->settings->cg_tol_reduction, "cg_tol_reduction");
  res.push_back(solver->settings->cg_tol_fraction, "cg_tol_fraction");
  res.push_back((int)solver->settings->cg_precond, "cg_precond");
  res.push_back((bool)solver->settings->check_dualgap, "check_dualgap");
  res.push_back(solver->settings->profiler_level, "profiler_level");

  // Deprecated aliases for backward compatibility
  res.push_back((bool)solver->settings->polishing, "polish");
  res.push_back((bool)solver->settings->warm_starting, "warm_start");

  return res;
}


// [[Rcpp::export]]
IntegerVector osqpGetDims(SEXP workPtr)
{
  auto solver = as<Rcpp::XPtr<OSQPSolver, Rcpp::PreserveStorage, mycleanup>>(workPtr);
  OSQPInt m_dim, n_dim;
  osqp_get_dimensions(solver, &m_dim, &n_dim);
  auto res = IntegerVector::create(_("n") = n_dim, _("m") = m_dim);
  return res;
}

// [[Rcpp::export]]
void osqpUpdate(SEXP workPtr,
    Rcpp::Nullable<NumericVector> q_new,
    Rcpp::Nullable<NumericVector> l_new,
    Rcpp::Nullable<NumericVector> u_new,
    Rcpp::Nullable<NumericVector> Px,
    Rcpp::Nullable<IntegerVector> Px_idx,
    Rcpp::Nullable<NumericVector> Ax,
    Rcpp::Nullable<IntegerVector> Ax_idx)
{
  auto solver = as<Rcpp::XPtr<OSQPSolver, Rcpp::PreserveStorage, mycleanup>>(workPtr);

  // Update problem vectors using consolidated function
  {
    const OSQPFloat* q_ptr = OSQP_NULL;
    const OSQPFloat* l_ptr = OSQP_NULL;
    const OSQPFloat* u_ptr = OSQP_NULL;

    // We need to keep vectors alive while osqp_update_data_vec uses them
    std::vector<OSQPFloat> q_vec, l_vec, u_vec;

    if (q_new.isNotNull()) {
      NumericVector q = q_new.get();
      q_vec.assign(q.begin(), q.end());
      q_ptr = q_vec.data();
    }
    if (l_new.isNotNull()) {
      NumericVector l = l_new.get();
      l_vec.assign(l.begin(), l.end());
      std::replace_if(l_vec.begin(), l_vec.end(), below_osqp_neg_inf, -OSQP_INFTY);
      std::replace_if(l_vec.begin(), l_vec.end(), above_osqp_inf, OSQP_INFTY);
      l_ptr = l_vec.data();
    }
    if (u_new.isNotNull()) {
      NumericVector u = u_new.get();
      u_vec.assign(u.begin(), u.end());
      std::replace_if(u_vec.begin(), u_vec.end(), below_osqp_neg_inf, -OSQP_INFTY);
      std::replace_if(u_vec.begin(), u_vec.end(), above_osqp_inf, OSQP_INFTY);
      u_ptr = u_vec.data();
    }

    if (q_ptr || l_ptr || u_ptr) {
      OSQPInt res = osqp_update_data_vec(solver, q_ptr, l_ptr, u_ptr);
      if (res != 0) stop("osqp_update_data_vec failed with error code %d", res);
    }
  }

  // Update problem matrices using consolidated function
  {
    const OSQPFloat* Px_ptr = OSQP_NULL;
    const OSQPInt*   Px_idx_ptr = OSQP_NULL;
    OSQPInt          P_new_n = 0;
    const OSQPFloat* Ax_ptr = OSQP_NULL;
    const OSQPInt*   Ax_idx_ptr = OSQP_NULL;
    OSQPInt          A_new_n = 0;

    std::vector<OSQPFloat> Px_vec, Ax_vec;
    std::vector<OSQPInt> Px_idx_vec, Ax_idx_vec;

    if (Px.isNotNull()) {
      NumericVector Px_r = Px.get();
      Px_vec.assign(Px_r.begin(), Px_r.end());
      Px_ptr = Px_vec.data();
      P_new_n = Px_vec.size();

      if (Px_idx.isNotNull()) {
        IntegerVector Px_idx_r = Px_idx.get();
        Px_idx_vec.assign(Px_idx_r.begin(), Px_idx_r.end());
        Px_idx_ptr = Px_idx_vec.data();
      }
    }

    if (Ax.isNotNull()) {
      NumericVector Ax_r = Ax.get();
      Ax_vec.assign(Ax_r.begin(), Ax_r.end());
      Ax_ptr = Ax_vec.data();
      A_new_n = Ax_vec.size();

      if (Ax_idx.isNotNull()) {
        IntegerVector Ax_idx_r = Ax_idx.get();
        Ax_idx_vec.assign(Ax_idx_r.begin(), Ax_idx_r.end());
        Ax_idx_ptr = Ax_idx_vec.data();
      }
    }

    if (Px_ptr || Ax_ptr) {
      OSQPInt res = osqp_update_data_mat(solver, Px_ptr, Px_idx_ptr, P_new_n, Ax_ptr, Ax_idx_ptr, A_new_n);
      if (res != 0) stop("osqp_update_data_mat failed with error code %d", res);
    }
  }
}

void extractMatrixData(const S4& mat, std::vector<OSQPInt>& iout, std::vector<OSQPInt>& pout, std::vector<OSQPFloat>& xout)
{
  IntegerVector i = mat.slot("i");
  IntegerVector p = mat.slot("p");
  NumericVector x = mat.slot("x");

  iout.resize(i.size());
  pout.resize(p.size());
  xout.resize(x.size());
  std::copy(i.begin(), i.end(), iout.begin());
  std::copy(p.begin(), p.end(), pout.begin());
  std::copy(x.begin(), x.end(), xout.begin());

  return;
}


void translateSettings(OSQPSettings* settings, const List& pars)
{
  CharacterVector nms(pars.names());
  for (int i = 0; i < pars.size(); i++)
  {
    if (Rf_isNull(nms[i]))
      continue;
    auto nm = as<std::string>(nms[i]);

    // Float settings
    if (nm == "rho")
      settings->rho = as<OSQPFloat>(pars[i]);
    else if (nm == "sigma")
      settings->sigma = as<OSQPFloat>(pars[i]);
    else if (nm == "eps_abs")
      settings->eps_abs = as<OSQPFloat>(pars[i]);
    else if (nm == "eps_rel")
      settings->eps_rel = as<OSQPFloat>(pars[i]);
    else if (nm == "eps_prim_inf")
      settings->eps_prim_inf = as<OSQPFloat>(pars[i]);
    else if (nm == "eps_dual_inf")
      settings->eps_dual_inf = as<OSQPFloat>(pars[i]);
    else if (nm == "alpha")
      settings->alpha = as<OSQPFloat>(pars[i]);
    else if (nm == "delta")
      settings->delta = as<OSQPFloat>(pars[i]);
    else if (nm == "adaptive_rho_fraction")
      settings->adaptive_rho_fraction = as<OSQPFloat>(pars[i]);
    else if (nm == "adaptive_rho_tolerance")
      settings->adaptive_rho_tolerance = as<OSQPFloat>(pars[i]);
    else if (nm == "time_limit")
      settings->time_limit = as<OSQPFloat>(pars[i]);
    else if (nm == "cg_tol_fraction")
      settings->cg_tol_fraction = as<OSQPFloat>(pars[i]);

    // Integer settings
    else if (nm == "linsys_solver")
      settings->linsys_solver = (osqp_linsys_solver_type)as<OSQPInt>(pars[i]);
    else if (nm == "polish_refine_iter")
      settings->polish_refine_iter = as<OSQPInt>(pars[i]);
    else if (nm == "check_termination")
      settings->check_termination = as<OSQPInt>(pars[i]);
    else if (nm == "scaling")
      settings->scaling = as<OSQPInt>(pars[i]);
    else if (nm == "max_iter")
      settings->max_iter = as<OSQPInt>(pars[i]);
    else if (nm == "adaptive_rho")
      settings->adaptive_rho = as<OSQPInt>(pars[i]);
    else if (nm == "adaptive_rho_interval")
      settings->adaptive_rho_interval = as<OSQPInt>(pars[i]);
    else if (nm == "cg_max_iter")
      settings->cg_max_iter = as<OSQPInt>(pars[i]);
    else if (nm == "cg_tol_reduction")
      settings->cg_tol_reduction = as<OSQPInt>(pars[i]);
    else if (nm == "cg_precond")
      settings->cg_precond = (osqp_precond_type)as<OSQPInt>(pars[i]);
    else if (nm == "profiler_level")
      settings->profiler_level = as<OSQPInt>(pars[i]);

    // Boolean settings (stored as OSQPInt)
    else if (nm == "polishing")
      settings->polishing = as<OSQPInt>(pars[i]);
    else if (nm == "verbose")
      settings->verbose = as<OSQPInt>(pars[i]);
    else if (nm == "scaled_termination")
      settings->scaled_termination = as<OSQPInt>(pars[i]);
    else if (nm == "warm_starting")
      settings->warm_starting = as<OSQPInt>(pars[i]);
    else if (nm == "rho_is_vec")
      settings->rho_is_vec = as<OSQPInt>(pars[i]);
    else if (nm == "check_dualgap")
      settings->check_dualgap = as<OSQPInt>(pars[i]);
  }

  return;
}

// [[Rcpp::export]]
void osqpUpdateSettings(SEXP workPtr, SEXP val, std::string nm)
{
  auto solver = as<Rcpp::XPtr<OSQPSolver, Rcpp::PreserveStorage, mycleanup>>(workPtr);

  // rho must be updated via osqp_update_rho, not osqp_update_settings
  if (nm == "rho") {
    OSQPInt res = osqp_update_rho(solver, as<OSQPFloat>(val));
    if (res != 0) stop("osqp_update_rho failed with error code %d", res);
    return;
  }

  // Copy-modify-update pattern: osqp_update_settings requires a complete OSQPSettings struct
  OSQPSettings new_settings = *(solver->settings);

  // Float settings
  if (nm == "eps_abs")
    new_settings.eps_abs = as<OSQPFloat>(val);
  else if (nm == "eps_rel")
    new_settings.eps_rel = as<OSQPFloat>(val);
  else if (nm == "eps_prim_inf")
    new_settings.eps_prim_inf = as<OSQPFloat>(val);
  else if (nm == "eps_dual_inf")
    new_settings.eps_dual_inf = as<OSQPFloat>(val);
  else if (nm == "alpha")
    new_settings.alpha = as<OSQPFloat>(val);
  else if (nm == "delta")
    new_settings.delta = as<OSQPFloat>(val);
  else if (nm == "time_limit")
    new_settings.time_limit = as<OSQPFloat>(val);
  else if (nm == "cg_tol_fraction")
    new_settings.cg_tol_fraction = as<OSQPFloat>(val);
  // Integer settings
  else if (nm == "max_iter")
    new_settings.max_iter = as<OSQPInt>(val);
  else if (nm == "check_termination")
    new_settings.check_termination = as<OSQPInt>(val);
  else if (nm == "polish_refine_iter")
    new_settings.polish_refine_iter = as<OSQPInt>(val);
  else if (nm == "cg_max_iter")
    new_settings.cg_max_iter = as<OSQPInt>(val);
  else if (nm == "cg_tol_reduction")
    new_settings.cg_tol_reduction = as<OSQPInt>(val);
  // Boolean settings
  else if (nm == "polishing")
    new_settings.polishing = as<OSQPInt>(val);
  else if (nm == "verbose")
    new_settings.verbose = as<OSQPInt>(val);
  else if (nm == "scaled_termination")
    new_settings.scaled_termination = as<OSQPInt>(val);
  else if (nm == "warm_starting")
    new_settings.warm_starting = as<OSQPInt>(val);
  else if (nm == "check_dualgap")
    new_settings.check_dualgap = as<OSQPInt>(val);
  else {
    Rcout << "Param " + nm + " cannot be updated live" << std::endl;
    return;
  }

  OSQPInt res = osqp_update_settings(solver, &new_settings);
  if (res != 0) stop("osqp_update_settings failed with error code %d", res);
}

// [[Rcpp::export]]
SEXP osqpGetData(SEXP workPtr, std::string nm)
{
  auto solver = as<Rcpp::XPtr<OSQPSolver, Rcpp::PreserveStorage, mycleanup>>(workPtr);

  // Access internal data via solver->work->data (private, but we include the headers)
  // Note: data is scaled after setup, so returned values are the scaled versions
  OSQPData* data = solver->work->data;

  if (nm == "P")
    return toDgCMat(data->P->csc);
  if (nm == "A")
    return toDgCMat(data->A->csc);

  if (nm == "q") {
    int n = data->n;
    NumericVector q(data->q->values, data->q->values + n);
    return q;
  }
  if (nm == "l") {
    int m = data->m;
    NumericVector l(data->l->values, data->l->values + m);
    return l;
  }
  if (nm == "u") {
    int m = data->m;
    NumericVector u(data->u->values, data->u->values + m);
    return u;
  }

  return R_NilValue;
}

S4 toDgCMat(OSQPCscMatrix* inmat)
{
  S4 m("dgCMatrix");

  int nnz = inmat->nzmax;
  int nr = inmat->m;
  int nc = inmat->n;

  NumericVector x(inmat->x, inmat->x + nnz);
  IntegerVector i(inmat->i, inmat->i + nnz);
  IntegerVector p(inmat->p, inmat->p + nc + 1);
  IntegerVector dim = IntegerVector::create(nr, nc);

  m.slot("i")   = i;
  m.slot("p")   = p;
  m.slot("x")   = x;
  m.slot("Dim") = dim;

  return m;
}


void mycleanup(OSQPSolver* x)
{
  osqp_cleanup(x);
}


// [[Rcpp::export]]
void osqpWarmStart(SEXP workPtr, Rcpp::Nullable<NumericVector> x, Rcpp::Nullable<NumericVector> y)
{
  auto solver = as<Rcpp::XPtr<OSQPSolver, Rcpp::PreserveStorage, mycleanup>>(workPtr);

  const OSQPFloat* x_ptr = OSQP_NULL;
  const OSQPFloat* y_ptr = OSQP_NULL;
  std::vector<OSQPFloat> x_vec, y_vec;

  if (x.isNotNull()) {
    NumericVector xr = x.get();
    x_vec.assign(xr.begin(), xr.end());
    x_ptr = x_vec.data();
  }
  if (y.isNotNull()) {
    NumericVector yr = y.get();
    y_vec.assign(yr.begin(), yr.end());
    y_ptr = y_vec.data();
  }

  if (x_ptr || y_ptr) {
    OSQPInt res = osqp_warm_start(solver, x_ptr, y_ptr);
    if (res != 0) stop("osqp_warm_start failed with error code %d", res);
  }
}

// [[Rcpp::export]]
void osqpColdStart(SEXP workPtr)
{
  auto solver = as<Rcpp::XPtr<OSQPSolver, Rcpp::PreserveStorage, mycleanup>>(workPtr);
  osqp_cold_start(solver);
}
