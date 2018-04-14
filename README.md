# R interface for OSQP

[![Build Status](https://travis-ci.org/goulart-paul/osqp-r)](https://travis-ci.org/goulart-paul/osqp-r.svg?branch=master)
[![Build Status](https://ci.appveyor.com/project/goulart-paul/osqp-r/)](https://ci.appveyor.com/project/goulart-paul/osqp-r/branch/master?svg=true)
  
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

From github:
```{r}
install.packages("remotes")
remotes::install_github("r-lib/remotes#103")
remotes::install_git("git://github.com/goulart-paul/osqp-r",submodules = TRUE)
```
