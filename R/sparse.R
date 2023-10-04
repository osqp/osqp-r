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
    Matrix::sparseMatrix(i = m$i, j = m$j, x = m$v, dims = dim(m))
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
  if (inherits(m, "dtCMatrix")) {
    m
  } else if (inherits(m, "matrix")) {
    Matrix::.m2sparse(m, "dtCMatrix")
  } else if(inherits(m, "simple_triplet_matrix")) { ## e.g. from package slam
    ind <- which(m$i <= m$j)
    Matrix::sparseMatrix(i = m$i[ind], j = m$j[ind], x = m$v[ind], dims = dim(m), triangular = TRUE)
  } else {
    ## Resort to brute force
    Matrix::triu(as(as(as(m, "CsparseMatrix"), "generalMatrix"), "dMatrix"))
  }
}
