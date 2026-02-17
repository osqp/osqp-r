# OSQP Solver object

OSQP Solver object

## Usage

``` r
osqp(P = NULL, q = NULL, A = NULL, l = NULL, u = NULL, pars = osqpSettings())
```

## Arguments

- P, A:

  sparse matrices of class dgCMatrix or coercible into such, with P
  positive semidefinite. (In the interest of efficiency, only the upper
  triangular part of P is used)

- q, l, u:

  Numeric vectors, with possibly infinite elements in l and u

- pars:

  list with optimization parameters, conveniently set with the function
  [`osqpSettings`](https://osqp.github.io/osqp-r/reference/osqpSettings.md).
  For `model@UpdateSettings(newPars)` only a subset of the settings can
  be updated once the problem has been initialized.

## Value

An S7 object of class "OSQP_Model" with computed properties that return
methods.

## Details

Allows one to solve a parametric problem with for example warm starts
between updates of the parameter, c.f. the examples. The object returned
by `osqp` contains several computed properties (accessed via `@`) which
can be used to either update/get details of the problem, modify the
optimization settings or attempt to solve the problem.

## Usage

    model = osqp(P=NULL, q=NULL, A=NULL, l=NULL, u=NULL, pars=osqpSettings())

    model@Solve()
    model@Update(q = NULL, l = NULL, u = NULL, Px = NULL, Px_idx = NULL, Ax = NULL, Ax_idx = NULL)
    model@GetParams()
    model@GetDims()
    model@UpdateSettings(newPars = list())

    model@GetData(element = c("P", "q", "A", "l", "u"))
    model@WarmStart(x=NULL, y=NULL)
    model@ColdStart()

    print(model)

## Method Arguments

- element:

  a string with the name of one of the matrices / vectors of the problem

- newPars:

  list with optimization parameters

## See also

[`solve_osqp`](https://osqp.github.io/osqp-r/reference/solve_osqp.md)

## Examples

``` r
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

settings <- osqpSettings(verbose = FALSE)

model <- osqp(P, q, A, l, u, settings)

# Solve
res <- model@Solve()

# Define new vector
q_new <- c(10., 20.)

# Update model and solve again
model@Update(q = q_new)
res <- model@Solve()
```
