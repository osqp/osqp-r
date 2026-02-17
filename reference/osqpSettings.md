# Settings for OSQP

For further details please consult the OSQP documentation:
<https://osqp.org/>

## Usage

``` r
osqpSettings(
  rho = 0.1,
  rho_is_vec = TRUE,
  sigma = 1e-06,
  max_iter = 4000L,
  eps_abs = 0.001,
  eps_rel = 0.001,
  eps_prim_inf = 1e-04,
  eps_dual_inf = 1e-04,
  alpha = 1.6,
  linsys_solver = c(OSQP_DIRECT_SOLVER = 1L),
  delta = 1e-06,
  polishing = FALSE,
  polish_refine_iter = 3L,
  verbose = TRUE,
  scaled_termination = FALSE,
  check_termination = 25L,
  check_dualgap = TRUE,
  warm_starting = TRUE,
  scaling = 10L,
  adaptive_rho = 1L,
  adaptive_rho_interval = 50L,
  adaptive_rho_tolerance = 5,
  adaptive_rho_fraction = 0.4,
  cg_max_iter = 20L,
  cg_tol_reduction = 10L,
  cg_tol_fraction = 0.15,
  cg_precond = c(OSQP_DIAGONAL_PRECONDITIONER = 1L),
  profiler_level = 0L,
  time_limit = 1e+10,
  polish = NULL,
  warm_start = NULL
)
```

## Arguments

- rho:

  ADMM step rho

- rho_is_vec:

  boolean, whether rho is treated as a vector (per-constraint) or scalar

- sigma:

  ADMM step sigma

- max_iter:

  maximum iterations

- eps_abs:

  absolute convergence tolerance

- eps_rel:

  relative convergence tolerance

- eps_prim_inf:

  primal infeasibility tolerance

- eps_dual_inf:

  dual infeasibility tolerance

- alpha:

  relaxation parameter

- linsys_solver:

  which linear systems solver to use, 1=OSQP_DIRECT_SOLVER (QDLDL),
  2=OSQP_INDIRECT_SOLVER

- delta:

  regularization parameter for polishing

- polishing:

  boolean, polish ADMM solution

- polish_refine_iter:

  iterative refinement steps in polishing

- verbose:

  boolean, write out progress

- scaled_termination:

  boolean, use scaled termination criteria

- check_termination:

  integer, check termination interval. If 0, termination checking is
  disabled

- check_dualgap:

  boolean, check duality gap termination criteria

- warm_starting:

  boolean, warm start

- scaling:

  heuristic data scaling iterations. If 0, scaling disabled

- adaptive_rho:

  integer, rho adaptation strategy: 0=disabled, 1=iterations, 2=time,
  3=KKT error

- adaptive_rho_interval:

  Number of iterations between rho adaptations rho. If 0, it is
  automatic

- adaptive_rho_tolerance:

  Tolerance X for adapting rho. The new rho has to be X times larger or
  1/X times smaller than the current one to trigger a new factorization

- adaptive_rho_fraction:

  Interval for adapting rho (fraction of the setup time)

- cg_max_iter:

  maximum number of CG iterations (indirect solver only)

- cg_tol_reduction:

  integer, number of consecutive zero CG iterations before the tolerance
  gets halved (indirect solver only)

- cg_tol_fraction:

  CG tolerance fraction (indirect solver only)

- cg_precond:

  preconditioner for CG method (indirect solver only): 0=none,
  1=diagonal (Jacobi)

- profiler_level:

  integer, level of detail for profiler annotations (0=off)

- time_limit:

  run time limit in seconds (1e10 effectively disables)

- polish:

  Deprecated. Use `polishing` instead.

- warm_start:

  Deprecated. Use `warm_starting` instead.
