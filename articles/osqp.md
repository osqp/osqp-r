# Solving Quadratic Programs with OSQP

## Introduction

The `osqp` package provides R bindings to the [OSQP](https://osqp.org/)
(Operator Splitting Quadratic Program) solver. OSQP solves convex
quadratic programs of the form $$\begin{array}{ll}
\text{minimize} & {\frac{1}{2}x^{T}Px + q^{T}x} \\
\text{subject to} & {l \leq Ax \leq u}
\end{array}$$ where $x \in \mathbf{R}^{n}$ is the optimization variable,
$P \in \mathbf{S}_{+}^{n}$ is a positive semidefinite matrix,
$q \in \mathbf{R}^{n}$, $A \in \mathbf{R}^{m \times n}$, and
$l,u \in \mathbf{R}^{m}$ (with elements possibly equal to $- \infty$ or
$+ \infty$).

The solver is based on the Alternating Direction Method of Multipliers
(ADMM). Key features include:

- **Efficient**: exploits structure in sparse QPs, requires only a
  single matrix factorization
- **Robust**: places no requirements on problem data (e.g., P need not
  be strictly positive definite); detects primal and dual infeasibility
- **Warm starting**: allows reusing previous solutions to speed up
  parametric solves
- **Polishing**: an optional step to obtain high-accuracy solutions

For algorithm details see [Stellato et al.
(2020)](https://doi.org/10.1007/s12532-020-00179-2).

## Migrating from osqp \< 1.0

Version 1.0.0 brings two major changes:

1.  **S7 classes replace R6.** The model object returned by
    [`osqp()`](https://osqp.github.io/osqp-r/reference/osqp.md) is now
    an S7 object. Methods are accessed with `@` instead of `$`:

    | Old (\< 1.0)                | New (\>= 1.0)               |
    |:----------------------------|:----------------------------|
    | `model$Solve()`             | `model@Solve()`             |
    | `model$Update(...)`         | `model@Update(...)`         |
    | `model$WarmStart(...)`      | `model@WarmStart(...)`      |
    | `model$ColdStart()`         | `model@ColdStart()`         |
    | `model$GetParams()`         | `model@GetParams()`         |
    | `model$GetDims()`           | `model@GetDims()`           |
    | `model$UpdateSettings(...)` | `model@UpdateSettings(...)` |
    | `model$GetData(...)`        | `model@GetData(...)`        |

    During the transition, the old `$` syntax still works but emits a
    deprecation warning. It will be removed in a future release.

2.  **OSQP C library upgraded from 0.6.3 to 1.0.0.** Some solver
    settings have been renamed:

    | Old (0.6.x)  | New (1.0)       |
    |:-------------|:----------------|
    | `warm_start` | `warm_starting` |
    | `polish`     | `polishing`     |

    The old names still work but emit a deprecation warning. They will
    be removed in a future release.

## Setup and Solve

We demonstrate the basic workflow with a small QP. First, load the
required packages.

``` r
library(osqp)
library(Matrix)
```

Define the problem data. The objective matrix `P` and constraint matrix
`A` should be sparse (from the `Matrix` package). Only the upper
triangular part of `P` is used.

``` r
P <- Matrix(c(4., 1.,
              1., 2.), 2, 2, sparse = TRUE)
q <- c(1., 1.)
A <- Matrix(c(1., 1., 0.,
              1., 0., 1.), 3, 2, sparse = TRUE)
l <- c(1., 0., 0.)
u <- c(1., 0.7, 0.7)
```

This corresponds to: $$\begin{array}{ll}
\text{minimize} & {\frac{1}{2}x^{T}\begin{bmatrix}
4 & 1 \\
1 & 2
\end{bmatrix}x + \begin{bmatrix}
1 \\
1
\end{bmatrix}^{T}x} \\
\text{subject to} & {\begin{bmatrix}
1 \\
0 \\
0
\end{bmatrix} \leq \begin{bmatrix}
1 & 1 \\
1 & 0 \\
0 & 1
\end{bmatrix}x \leq \begin{bmatrix}
1 \\
0.7 \\
0.7
\end{bmatrix}}
\end{array}$$

Create the model and solve.

``` r
model <- osqp(P, q, A, l, u, pars = osqpSettings(verbose = FALSE))
result <- model@Solve()
```

The result is a list containing the primal solution `x`, dual solution
`y`, and solver information in `info`.

``` r
result$x
#> [1] 0.3013757 0.6983957
result$y
#> [1] -2.904627e+00 -1.385901e-18  2.055120e-01
result$info$status
#> [1] "solved"
result$info$obj_val
#> [1] 1.879662
```

### Convenience Function

For one-shot solves where you do not need to update the problem
afterward,
[`solve_osqp()`](https://osqp.github.io/osqp-r/reference/solve_osqp.md)
provides a simpler interface.

``` r
result <- solve_osqp(P, q, A, l, u, pars = osqpSettings(verbose = FALSE))
result$x
#> [1] 0.3013757 0.6983957
```

## Updating Problem Data

A major advantage of OSQP is the ability to update problem vectors and
matrices *without* re-running the full setup. This avoids re-factorizing
the KKT system and is much faster than creating a new model from
scratch.

### Updating Vectors

The vectors $q$, $l$, and $u$ can be updated freely.

``` r
q_new <- c(2., 3.)
l_new <- c(2., -1., -1.)
u_new <- c(2., 2.5, 2.5)
model@Update(q = q_new, l = l_new, u = u_new)
result <- model@Solve()
result$x
#> [1] 0.7499527 1.2499010
result$info$status
#> [1] "solved"
```

### Updating Matrices

The nonzero values of $P$ and $A$ can be updated as long as the
**sparsity pattern remains the same**. Pass the new nonzero values via
`Px` and `Ax`.

``` r
# Create new matrices with same sparsity pattern
P_new <- Matrix(c(5., 1.5,
                  1.5, 1.), 2, 2, sparse = TRUE)
A_new <- Matrix(c(1.2, 1.5, 0.,
                  1.1, 0., 0.8), 3, 2, sparse = TRUE)

# Update only the nonzero values
model@Update(Px = triu(P_new)@x, Ax = A_new@x)
result <- model@Solve()
result$x
#> [1] 0.1813419 1.6204746
```

Note that for `P` you must pass only the **upper triangular** nonzero
values (using
[`triu()`](https://rdrr.io/pkg/Matrix/man/band-methods.html)), since
OSQP stores $P$ in upper triangular form.

## Warm Starting

When solving a sequence of related problems (e.g., in a parameter sweep
or iterative algorithm), warm starting can dramatically reduce the
number of ADMM iterations by initializing the solver with a good
starting point.

``` r
# Solve original problem
model <- osqp(P, q, A, l, u, pars = osqpSettings(verbose = FALSE,
                                                   warm_starting = TRUE))
res1 <- model@Solve()
cold_iters <- res1$info$iter

# Warm start with previous solution
model@WarmStart(x = res1$x, y = res1$y)
res2 <- model@Solve()
warm_iters <- res2$info$iter

cat("Cold start iterations:", cold_iters, "\n")
#> Cold start iterations: 25
cat("Warm start iterations:", warm_iters, "\n")
#> Warm start iterations: 25
```

You can warm start with just the primal variables (`x`), just the dual
variables (`y`), or both. To reset the iterate and start from scratch,
use `ColdStart()`.

``` r
model@ColdStart()
res3 <- model@Solve()
cat("After cold start:", res3$info$iter, "iterations\n")
#> After cold start: 25 iterations
```

## Solver Settings

Solver behavior is controlled through
[`osqpSettings()`](https://osqp.github.io/osqp-r/reference/osqpSettings.md).
Settings are passed at model creation time, and a subset can be updated
afterward.

``` r
# Tighter tolerances and solution polishing
settings <- osqpSettings(
  eps_abs   = 1e-06,
  eps_rel   = 1e-06,
  polishing = TRUE,
  verbose   = FALSE
)
model <- osqp(P, q, A, l, u, pars = settings)
result <- model@Solve()
result$info$obj_val
#> [1] 1.88
```

### Key Settings

| Setting             | Default | Description                                         |
|:--------------------|:-------:|:----------------------------------------------------|
| `eps_abs`           |  1e-3   | Absolute convergence tolerance                      |
| `eps_rel`           |  1e-3   | Relative convergence tolerance                      |
| `max_iter`          |  4000   | Maximum ADMM iterations                             |
| `polishing`         |  FALSE  | Polish the ADMM solution for higher accuracy        |
| `warm_starting`     |  TRUE   | Enable warm starting                                |
| `verbose`           |  TRUE   | Print solver output                                 |
| `scaling`           |   10    | Data scaling iterations (0 to disable)              |
| `adaptive_rho`      |    1    | Adaptive ADMM step size (0=off, 1=iterations)       |
| `time_limit`        |  1e10   | Time limit in seconds                               |
| `check_termination` |   25    | Check termination every N iterations (0 to disable) |

See
[`?osqpSettings`](https://osqp.github.io/osqp-r/reference/osqpSettings.md)
for the full list.

### Updating Settings

Some settings can be updated after the model is created.

``` r
model@UpdateSettings(osqpSettings(max_iter = 2000L))
pars <- model@GetParams()
pars$max_iter
#> [1] 2000
```

### Retrieving Problem Dimensions

``` r
dims <- model@GetDims()
cat("Variables:", dims[["n"]], " Constraints:", dims[["m"]], "\n")
#> Variables: 2  Constraints: 3
```

## Result Object

The result from `Solve()` contains:

| Field           | Description                                      |
|:----------------|:-------------------------------------------------|
| `x`             | Primal solution ($n$-vector)                     |
| `y`             | Dual solution ($m$-vector)                       |
| `prim_inf_cert` | Primal infeasibility certificate (if applicable) |
| `dual_inf_cert` | Dual infeasibility certificate (if applicable)   |
| `info`          | List of solver information                       |

The `info` list includes:

| Field        | Description                      |
|:-------------|:---------------------------------|
| `iter`       | Number of ADMM iterations        |
| `status`     | Status string (e.g., `"solved"`) |
| `status_val` | Status integer code              |
| `obj_val`    | Optimal objective value          |
| `prim_res`   | Primal residual                  |
| `dual_res`   | Dual residual                    |
| `setup_time` | Time for setup (seconds)         |
| `solve_time` | Time for solve (seconds)         |
| `run_time`   | Total time (seconds)             |

### Status Values

| Status                     | Constant                  | Value |
|:---------------------------|:--------------------------|:-----:|
| solved                     | `OSQP_SOLVED`             |   1   |
| solved inaccurate          | `OSQP_SOLVED_INACCURATE`  |   2   |
| primal infeasible          | `OSQP_PRIMAL_INFEASIBLE`  |  -3   |
| dual infeasible            | `OSQP_DUAL_INFEASIBLE`    |  -4   |
| maximum iterations reached | `OSQP_MAX_ITER_REACHED`   |  -2   |
| run time limit reached     | `OSQP_TIME_LIMIT_REACHED` |  -6   |

## Example: Lasso Regression

[Lasso](https://en.wikipedia.org/wiki/Lasso_(statistics)) performs
sparse linear regression by adding an $\ell_{1}$ penalty:
$$\text{minimize}\quad\frac{1}{2} \parallel A_{d}x - b \parallel_{2}^{2} + \gamma \parallel x \parallel_{1}$$

This can be reformulated as a QP by introducing auxiliary variables $y$
and $t$: $$\begin{array}{ll}
\text{minimize} & {\frac{1}{2}y^{T}y + \gamma\mathbf{1}^{T}t} \\
\text{subject to} & {y = A_{d}x - b} \\
 & {- t \leq x \leq t}
\end{array}$$

Because the regularization parameter $\gamma$ appears only in the linear
cost $q$, we can update it without re-factorizing the KKT system. Warm
starting makes sweeping over $\gamma$ very efficient.

``` r
set.seed(1)
n <- 10   # number of features
m <- 100  # number of observations

# Generate random data
Ad <- Matrix(rnorm(m * n), m, n, sparse = FALSE)
Ad <- as(Ad, "CsparseMatrix")
x_true <- rnorm(n) * (runif(n) > 0.8) / sqrt(n)
b <- as.numeric(Ad %*% x_true + 0.5 * rnorm(m))

# Auxiliary matrices
In <- Diagonal(n)
Im <- Diagonal(m)
On <- Matrix(0, n, n, sparse = TRUE)
Onm <- Matrix(0, n, m, sparse = TRUE)

# QP formulation: variables are (x, y, t)
P <- bdiag(On, Im, On)
q <- rep(0, 2 * n + m)
A <- rbind(
  cbind(Ad, -Im, Matrix(0, m, n, sparse = TRUE)),
  cbind(In, Onm, -In),
  cbind(In, Onm,  In)
)
l <- c(b, rep(-Inf, n), rep(0, n))
u <- c(b, rep(0, n), rep(Inf, n))

# Setup model once
model <- osqp(P, q, A, l, u,
              pars = osqpSettings(warm_starting = TRUE, verbose = FALSE))

# Solve for a range of gamma values
gammas <- seq(0.1, 2.0, length.out = 10)
results <- data.frame(gamma = gammas, nnz = integer(10), obj = numeric(10))

for (i in seq_along(gammas)) {
  gamma <- gammas[i]
  q_new <- c(rep(0, n + m), rep(gamma, n))
  model@Update(q = q_new)
  res <- model@Solve()
  x_hat <- res$x[1:n]
  results$nnz[i] <- sum(abs(x_hat) > 1e-4)
  results$obj[i] <- res$info$obj_val
}

results
#>        gamma nnz      obj
#> 1  0.1000000  10 11.37229
#> 2  0.3111111  10 11.50198
#> 3  0.5222222  10 11.63016
#> 4  0.7333333  10 11.75222
#> 5  0.9444444   9 11.86987
#> 6  1.1555556  10 11.98151
#> 7  1.3666667  10 12.09577
#> 8  1.5777778  10 12.19795
#> 9  1.7888889   9 12.30656
#> 10 2.0000000   9 12.41062
```

## Example: Portfolio Optimization

Mean-variance portfolio optimization allocates assets to maximize
risk-adjusted return: $$\begin{array}{ll}
\text{maximize} & {\mu^{T}x - \gamma\, x^{T}\Sigma x} \\
\text{subject to} & {\mathbf{1}^{T}x = 1} \\
 & {x \geq 0}
\end{array}$$ where $x \in \mathbf{R}^{n}$ is the portfolio, $\mu$ the
expected returns, $\gamma > 0$ the risk aversion, and $\Sigma$ the
covariance matrix. Using a factor model $\Sigma = FF^{T} + D$ (with
$F \in \mathbf{R}^{n \times k}$, $D$ diagonal), we can reformulate this
as: $$\begin{array}{ll}
\text{minimize} & {\frac{1}{2}x^{T}Dx + \frac{1}{2}y^{T}y - \frac{1}{2\gamma}\mu^{T}x} \\
\text{subject to} & {y = F^{T}x} \\
 & {\mathbf{1}^{T}x = 1} \\
 & {x \geq 0}
\end{array}$$

``` r
set.seed(42)
n <- 20   # assets
k <- 5    # factors
gamma <- 1

# Random factor model
F_mat <- Matrix(rnorm(n * k) * (runif(n * k) > 0.3), n, k, sparse = TRUE)
D <- Diagonal(n, runif(n) * sqrt(k))
mu <- rnorm(n)

# QP data: variables are (x, y)
P <- bdiag(D, Diagonal(k))
q <- c(-mu / (2 * gamma), rep(0, k))
A <- rbind(
  cbind(t(F_mat), -Diagonal(k)),                          # y = F'x
  cbind(Matrix(1, 1, n, sparse = TRUE), Matrix(0, 1, k, sparse = TRUE)),  # 1'x = 1
  cbind(Diagonal(n), Matrix(0, n, k, sparse = TRUE))       # x >= 0
)
l <- c(rep(0, k), 1, rep(0, n))
u <- c(rep(0, k), 1, rep(1, n))

result <- solve_osqp(P, q, A, l, u,
                     pars = osqpSettings(verbose = FALSE, polishing = TRUE))
cat("Status:", result$info$status, "\n")
#> Status: solved

# Portfolio weights (first n variables)
weights <- result$x[1:n]
cat("Invested assets:", sum(weights > 1e-4), "of", n, "\n")
#> Invested assets: 8 of 20
cat("Sum of weights:", round(sum(weights), 6), "\n")
#> Sum of weights: 1
```

## Sparse Matrix Input

OSQP requires sparse matrices in compressed sparse column (CSC) format.
The package accepts several input types and coerces them automatically:

- **`dgCMatrix`** (from `Matrix`) — used directly
- **Dense `matrix`** — converted to sparse
- **`simple_triplet_matrix`** (from `slam`) — converted to sparse

The objective matrix `P` is stored as upper triangular internally; the
package handles this conversion for you.

``` r
# Dense matrix input works fine
P_dense <- matrix(c(4, 1, 1, 2), 2, 2)
q <- c(1, 1)
A_dense <- matrix(c(1, 1, 0, 1, 0, 1), 3, 2)
l <- c(1, 0, 0)
u <- c(1, 0.7, 0.7)

result <- solve_osqp(P_dense, q, A_dense, l, u,
                     pars = osqpSettings(verbose = FALSE))
result$x
#> [1] 0.3013757 0.6983957
```

## References

If you use OSQP in your work, please cite the following references.

- B. Stellato, G. Banjac, P. Goulart, A. Bemporad, and S. Boyd, “OSQP:
  An Operator Splitting Solver for Quadratic Programs,” *Mathematical
  Programming Computation*, 2020.
  [doi:10.1007/s12532-020-00179-2](https://doi.org/10.1007/s12532-020-00179-2)

- G. Banjac, P. Goulart, B. Stellato, and S. Boyd, “Infeasibility
  Detection in the Alternating Direction Method of Multipliers for
  Convex Optimization,” *Journal of Optimization Theory and
  Applications*, 2019.

Please also cite the R package: run `citation("osqp")` for details.

- B. Stellato, G. Banjac, P. Goulart, S. Boyd, V. Bansal, and B.
  Narasimhan, *osqp: Quadratic Programming Solver using the ‘OSQP’
  Library*, R package, <https://osqp.org>.
