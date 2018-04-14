# R interface for OSQP

[![Build Status](https://travis-ci.org/oxfordcontrol/osqp-r.svg?branch=master)](https://travis-ci.org/oxfordcontrol/osqp-r)
[![Build status](https://ci.appveyor.com/api/projects/status/bx1navxa474nhlpd/branch/master?svg=true)](https://ci.appveyor.com/project/goulart-paul/osqp-r/branch/master)
[![](https://www.r-pkg.org/badges/version/rosqp)](https://www.r-pkg.org/pkg/rosqp)
[![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/rosqp)](https://www.r-pkg.org/pkg/rosqp)

Provides R-bindings to the [OSQP-solver](http://osqp.readthedocs.io/)

OSQP is a sparse Quadratic Programming Solver suitable for large problems.
It solves problems of the form:
```
minimize        0.5 x' P x + q' x

subject to      l <= A x <= u
```

with P Positive semidefinite

### OSQP License
The source code of OSQP is provided in the package and it is licensed under Apache 2.0.
Copyright (c) 2017 Bartolomeo Stellato, Goran Banjac, Paul Goulart, Stephen Boyd

This product includes software developed at Stanford University and at the University of Oxford.

Please see the file LICENSE and inst/COPYRIGHT for more details.

### Example
```{r}
library(rosqp)
## example, adapted from ?quadprog::solve.QP
Dmat       <- diag(3)
dvec       <- c(0,-5,0)
Amat       <- matrix(c(-4, 2, 0, -3, 1, -2, 0, 0, 1),3,3)
bvec       <- c(-8,2,0)
res = solve_osqp(Dmat, dvec, Amat, bvec)
print(res$x)
```

### Installation

A pre-compiled version of the OSQP-R interface can be installed directly from within R.   Note that this will install the OSQP interface from the current CRAN repository, which may not be the most up-to-date version:

```{r}
install.packages("rosqp")
```

If you would like to use the most recent version of OSQP-R and have access to git on your machine, then you can do the following from within a terminal:

```
git clone --recursive https://github.com/oxfordcontrol/osqp-r.git
cd osqp-r
R CMD install .
```

If you would like to install the latest version and do not have access to git on your machine, then you can do the following from within R:

```{r}
install.packages("remotes")
remotes::install_github("r-lib/remotes#103")
remotes::install_git("git://github.com/OxfordControl/osqp-r",submodules = TRUE)
```

Note that the second line above is necessary because the "remotes" package in R does not currently support recursive cloning of git submodules.
