# Sparse Quadratic Programming Solver

Solves \$\$arg\min_x 0.5 x'P x + q'x\$\$ s.t. \$\$l_i \< (A x)\_i \<
u_i\$\$ for real matrices P (nxn, positive semidefinite) and A (mxn)
with m number of constraints

## Usage

``` r
solve_osqp(
  P = NULL,
  q = NULL,
  A = NULL,
  l = NULL,
  u = NULL,
  pars = osqpSettings()
)
```

## Arguments

- P, A:

  sparse matrices of class dgCMatrix or coercible into such, with P
  positive semidefinite. Only the upper triangular part of P will be
  used.

- q, l, u:

  Numeric vectors, with possibly infinite elements in l and u

- pars:

  list with optimization parameters, conveniently set with the function
  `osqpSettings`

## Value

A list with elements x (the primal solution), y (the dual solution),
prim_inf_cert, dual_inf_cert, and info.

## References

Stellato, B., Banjac, G., Goulart, P., Bemporad, A., Boyd and S. (2018).
“OSQP: An Operator Splitting Solver for Quadratic Programs.” *ArXiv
e-prints*. 1711.08013.

## See also

[`osqp`](https://osqp.github.io/osqp-r/reference/osqp.md). The
underlying OSQP documentation: <https://osqp.org/>

## Examples

``` r
library(osqp)
## example, adapted from OSQP documentation
library(Matrix)

P <- Matrix(c(11., 0.,
              0., 0.), 2, 2, sparse = TRUE)
q <- c(3., 4.)
A <- Matrix(c(-1., 0., -1., 2., 3.,
              0., -1., -3., 5., 4.)
              , 5, 2, sparse = TRUE)
u <- c(0., 0., -15., 100., 80)
l <- rep_len(-Inf, 5)

settings <- osqpSettings(verbose = TRUE)

# Solve with OSQP
res <- solve_osqp(P, q, A, l, u, settings)
#> -----------------------------------------------------------------
#>            OSQP v1.0.0  -  Operator Splitting QP Solver
#>               (c) The OSQP Developer Team
#> -----------------------------------------------------------------
#> problem:  variables n = 2, constraints m = 5
#>           nnz(P) + nnz(A) = 9
#> settings: algebra = Built-in,
#>           OSQPInt = 4 bytes, OSQPFloat = 8 bytes,
#>           linear system solver = QDLDL v0.1.8,
#>           eps_abs = 1.0e-03, eps_rel = 1.0e-03,
#>           eps_prim_inf = 1.0e-04, eps_dual_inf = 1.0e-04,
#>           rho = 1.00e-01 (adaptive: 50 iterations),
#>           sigma = 1.00e-06, alpha = 1.60, max_iter = 4000
#>           check_termination: on (interval 25, duality gap: on),
#>           time_limit: 1.00e+10 sec,
#>           scaling: on (10 iterations), scaled_termination: off
#>           warm starting: on, polishing: off, 
#> iter   objective    prim res   dual res   gap        rel kkt    rho         time
#>    1  -7.2491e+00   2.04e+01   3.59e+00  -3.73e+01   2.04e+01   1.00e-01    3.28e-05s
#>  100   2.0002e+01   1.76e-03   1.47e-04   1.61e-03   1.76e-03   1.00e-01    6.71e-05s
#> 
#> status:               solved
#> number of iterations: 100
#> optimal objective:    20.0023
#> dual objective:       20.0007
#> duality gap:          1.6060e-03
#> primal-dual integral: 4.1060e+01
#> run time:             8.25e-05s
#> optimal rho estimate: 1.78e-01
#> 
res$x
#> [1] -3.125718e-08  5.000585e+00
```
