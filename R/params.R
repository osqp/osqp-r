
#' Settings for OSQP
#'
#' For further details please consult the OSQP documentation:
#' \url{https://osqp.org/}
#' @param rho ADMM step rho
#' @param rho_is_vec boolean, whether rho is treated as a vector (per-constraint) or scalar
#' @param sigma ADMM step sigma
#' @param max_iter maximum iterations
#' @param eps_rel relative convergence tolerance
#' @param eps_abs absolute convergence tolerance
#' @param eps_prim_inf primal infeasibility tolerance
#' @param eps_dual_inf dual infeasibility tolerance
#' @param alpha relaxation parameter
#' @param linsys_solver which linear systems solver to use, 1=OSQP_DIRECT_SOLVER (QDLDL), 2=OSQP_INDIRECT_SOLVER
#' @param delta regularization parameter for polishing
#' @param polishing boolean, polish ADMM solution
#' @param polish_refine_iter iterative refinement steps in polishing
#' @param verbose boolean, write out progress
#' @param scaled_termination boolean, use scaled termination criteria
#' @param check_termination integer, check termination interval. If 0, termination checking is disabled
#' @param check_dualgap boolean, check duality gap termination criteria
#' @param warm_starting boolean, warm start
#' @param scaling heuristic data scaling iterations. If 0, scaling disabled
#' @param adaptive_rho integer, rho adaptation strategy: 0=disabled, 1=iterations, 2=time, 3=KKT error
#' @param adaptive_rho_interval Number of iterations between rho adaptations rho. If 0, it is automatic
#' @param adaptive_rho_tolerance Tolerance X for adapting rho. The new rho has to be X times larger or 1/X times smaller than the current one to trigger a new factorization
#' @param adaptive_rho_fraction Interval for adapting rho (fraction of the setup time)
#' @param cg_max_iter maximum number of CG iterations (indirect solver only)
#' @param cg_tol_reduction integer, number of consecutive zero CG iterations before the tolerance gets halved (indirect solver only)
#' @param cg_tol_fraction CG tolerance fraction (indirect solver only)
#' @param cg_precond preconditioner for CG method (indirect solver only): 0=none, 1=diagonal (Jacobi)
#' @param profiler_level integer, level of detail for profiler annotations (0=off)
#' @param time_limit run time limit in seconds (1e10 effectively disables)
#' @param polish Deprecated. Use \code{polishing} instead.
#' @param warm_start Deprecated. Use \code{warm_starting} instead.
#' @export
osqpSettings = function(rho = 0.1, rho_is_vec = TRUE, sigma = 1e-06, max_iter = 4000L, eps_abs = 0.001,
                        eps_rel = 0.001, eps_prim_inf = 1e-04, eps_dual_inf = 1e-04,
                        alpha = 1.6, linsys_solver = c(OSQP_DIRECT_SOLVER=1L),
                        delta = 1e-06, polishing = FALSE, polish_refine_iter = 3L, verbose = TRUE,
                        scaled_termination = FALSE, check_termination = 25L, check_dualgap = TRUE,
                        warm_starting = TRUE,
                        scaling = 10L, adaptive_rho = 1L, adaptive_rho_interval = 50L,
                        adaptive_rho_tolerance = 5, adaptive_rho_fraction = 0.4,
                        cg_max_iter = 20L, cg_tol_reduction = 10L, cg_tol_fraction = 0.15,
                        cg_precond = c(OSQP_DIAGONAL_PRECONDITIONER=1L),
                        profiler_level = 0L,
                        time_limit = 1e10,
                        polish = NULL,
                        warm_start = NULL) {

  ## Deprecated aliases: polish -> polishing, warm_start -> warm_starting
  if (!is.null(polish)) {
    warning(cli::format_warning(
      c("!" = "The {.arg polish} argument is deprecated.",
        "i" = "Use {.arg polishing} instead.")),
      call. = FALSE
    )
    polishing <- polish
  }
  if (!is.null(warm_start)) {
    warning(cli::format_warning(
      c("!" = "The {.arg warm_start} argument is deprecated.",
        "i" = "Use {.arg warm_starting} instead.")),
      call. = FALSE
    )
    warm_starting <- warm_start
  }

  given_args  <- as.list(environment())      ## all params with current values
  call_arg_names  <- names(match.call()[-1]) ## lose the function name at index 1
  default_args  <- formals()                 ## this is the default list of arg values

  ## Map deprecated names to their replacements in the args list
  if ("polish" %in% call_arg_names) {
    call_arg_names <- c(setdiff(call_arg_names, "polish"), "polishing")
  }
  if ("warm_start" %in% call_arg_names) {
    call_arg_names <- c(setdiff(call_arg_names, "warm_start"), "warm_starting")
  }

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
