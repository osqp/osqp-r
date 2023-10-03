## Created by @naras to address coercion considerations with Matrix package
##
## Inspired by Rsymphony package from
## Kurt Hornik and Stefan Theussl and Reinhard Harter
##

#' Ensure that a matrix is column-sparse format (class `dgCMatrix`)
#' @param m a C sparse matrix or coercible object such as a matrix or simple triplet matrix
#' @return a sparse matrix of class `dgCMatrix`
#' @importFrom Matrix sparseMatrix
#' @importFrom methods as
#' @noRd
ensure_dgc_matrix <- function(m) {
  if (inherits(m, "dgCMatrix")) {
    m
  } else if (inherits(m, "matrix")) {
    Matrix::.m2sparse(m, "dgCMatrix")
  } else if(inherits(m, "simple_triplet_matrix")) { ## e.g. from package slam
    ## The matrix method assumes that indices for non-zero entries are
    ## in row-major order, but the simple_triplet_matrix() constructor
    ## currently does not canonicalize accordingly ...
    ind <- order(m$j, m$i)
    Matrix::sparseMatrix(p = c(0L, cumsum(tabulate(m$j[ind], m$ncol))),
                         i = m$i[ind] - 1L,
                         values = m$v[ind], dims = dim(m), index1 = FALSE)
  } else {
    ## Resort to brute force
    as(as(as(m, "CsparseMatrix"), "generalMatrix"), "dMatrix")
  }
}

#' Ensure that a matrix is column-sparse format upper triangular (class `dtCMatrix`)
#' @param m a C sparse upper triangular matrix or coercible object such as a matrix or simple triplet matrix
#' @return a sparse matrix of class `dgCMatrix`
#' @importFrom Matrix sparseMatrix triu
#' @importFrom methods as
#' @noRd
ensure_dtc_matrix <- function(m) {
  if (inherits(m, "dgCMatrix")) {
    m
  } else if (inherits(m, "matrix")) {
    Matrix::.m2sparse(m, "dtCMatrix")
  } else if(inherits(m, "simple_triplet_matrix")) { ## e.g. from package slam
    ind <- which(m$i <= m$j)
    x <- list(i = m$i[ind] + 1L, j = m$j[ind] + 1L) ##make it 1-based
    values <- m$v[ind]
    ind  <- order(x$j, x$i)  ## may not be needed
    Matrix::sparseMatrix(p = c(0L, cumsum(tabulate(x$j[ind], m$ncol))),
                         i = x$i[ind] - 1L,
                         values = values,
                         dims = dim(x), index1 = FALSE,
                         triangular = TRUE)
  } else {
    ## Resort to brute force
    Matrix::triu(as(as(as(m, "CsparseMatrix"), "generalMatrix"), "dMatrix"))
  }
}
