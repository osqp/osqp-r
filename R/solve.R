#' Sparse Quadratic Programming Solver
#' 
#' Solves \deqn{arg\min_x 0.5 x'P x + q'x}{argmin_x 0.5 x'P x + q'x}
#' s.t. \deqn{l_i < (A x)_i < u_i}{li < (A x)i < ui}
#' for real matrices P (nxn, positive semidefinite) and A (mxn) with m number of constraints
#' @param P,A sparse matrices of class dgCMatrix or coercible into such, with P positive semidefinite.
#' @param q,l,u Numeric vectors, with possibly infinite elements in l and u
#' @param pars list with optimization parameters, conveniently set with the function \code{osqpSettings} 
#' @return A list with elements x (the primal solution), y (the dual solution), prim_inf_cert, 
#' dual_inf_cert, and info.
#' @seealso \code{\link{osqp}}. The underlying OSQP documentation: \url{http://osqp.readthedocs.io/}
#' @examples
#' library(rosqp)
#' ## example, adapted from ?quadprog::solve.QP
#' Dmat       <- diag(3)
#' dvec       <- c(0,-5,0)
#' Amat       <- matrix(c(-4, 2, 0, -3, 1, -2, 0, 0, 1),3,3)
#' bvec       <- c(-8,2,0)
#' res = solve_osqp(Dmat, dvec, Amat, bvec)
#' print(res$x)
#' 
#' 
#' 
#' @references{
#' Stellato, B., Banjac, G., Goulart, P., Bemporad, A., Boyd and S. (2017).
#' ``OSQP: An Operator Splitting Solver for Quadratic Programs.''
#' \emph{ArXiv e-prints}.
#' 1711.08013.}
#' @export
solve_osqp <- function(P=NULL, q=NULL, A=NULL, l=NULL, u=NULL, pars = osqpSettings()) {
  
  model = osqp(P, q, A, l, u, pars)
  model$Solve()
}

