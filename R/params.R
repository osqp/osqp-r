
#' Settings for OSQP
#'
#' For further details please consult the OSQP documentation:
#' \url{https://osqp.org/}
#' @param rho ADMM step rho
#' @param sigma ADMM step sigma
#' @param max_iter maximum iterations
#' @param eps_rel relative convergence tolerance
#' @param eps_abs absolute convergence tolerance
#' @param eps_prim_inf primal infeasibility tolerance
#' @param eps_dual_inf dual infeasibility tolerance
#' @param alpha relaxation parameter
#' @param linsys_solver which linear systems solver to use, 0=QDLDL, 1=MKL Pardiso
#' @param delta regularization parameter for polish
#' @param polish boolean, polish ADMM solution
#' @param polish_refine_iter iterative refinement steps in polish
#' @param verbose boolean, write out progress
#' @param scaled_termination boolean, use scaled termination criteria
#' @param check_termination integer, check termination interval. If 0, termination checking is disabled
#' @param warm_start boolean, warm start
#' @param scaling heuristic data scaling iterations. If 0, scaling disabled
#' @param adaptive_rho cboolean, is rho step size adaptive?
#' @param adaptive_rho_interval Number of iterations between rho adaptations rho. If 0, it is automatic
#' @param adaptive_rho_tolerance Tolerance X for adapting rho. The new rho has to be X times larger or 1/X times smaller than the current one to trigger a new factorization
#' @param adaptive_rho_fraction Interval for adapting rho (fraction of the setup time)
#' @param time_limit run time limit with 0 indicating no limit
#' @export
osqpSettings = function(rho = 0.1, sigma = 1e-06, max_iter = 4000L, eps_abs = 0.001,
                        eps_rel = 0.001, eps_prim_inf = 1e-04, eps_dual_inf = 1e-04,
                        alpha = 1.6, linsys_solver = c(QDLDL_SOLVER=0L),
                        delta = 1e-06, polish = FALSE, polish_refine_iter = 3L, verbose = TRUE,
                        scaled_termination = FALSE, check_termination = 25L, warm_start = TRUE,
                        scaling = 10L, adaptive_rho = 1L, adaptive_rho_interval = 0L,
                        adaptive_rho_tolerance = 5, adaptive_rho_fraction = 0.4, time_limit = 0.0) {
  given_args  <- as.list(environment())      ## all params with current values
  call_arg_names  <- names(match.call()[-1]) ## lose the function name at index 1
  default_args  <- formals()                 ## this is the default list of arg values
  given_args  <- given_args[call_arg_names]  ## restrict to specified args

  for (name in call_arg_names) {
    given <- given_args[[name]]
    if (length(given) != 1 || is.na(given)) {
      given_args[[name]] <- default_args[[name]] ## force default
    } else {
      storage.mode(given_args[[name]]) <- storage.mode(eval(default_args[[name]])) #eval default arg
    }
  }
  given_args
}
