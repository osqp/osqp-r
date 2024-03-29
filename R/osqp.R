
#' OSQP Solver object
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom Matrix triu
#' @importFrom methods as
#' @importFrom R6 R6Class
#' @param P,A sparse matrices of class dgCMatrix or coercible into such, with P positive semidefinite. (In the interest of efficiency, only the upper triangular part of P is used)
#' @param q,l,u Numeric vectors, with possibly infinite elements in l and u
#' @param pars list with optimization parameters, conveniently set with the function
#' \code{\link{osqpSettings}}. For \code{osqpObject$UpdateSettings(newPars)} only a subset of the settings
#' can be updated once the problem has been initialized.
#' @seealso \code{\link{solve_osqp}}
#' @section Usage:
#' \preformatted{model = osqp(P=NULL, q=NULL, A=NULL, l=NULL, u=NULL, pars=osqpSettings())
#'
#' model$Solve()
#' model$Update(q = NULL, l = NULL, u = NULL, Px = NULL, Px_idx = NULL, Ax = NULL, Ax_idx = NULL)
#' model$GetParams()
#' model$GetDims()
#' model$UpdateSettings(newPars = list())
#'
#' model$GetData(element = c("P", "q", "A", "l", "u"))
#' model$WarmStart(x=NULL, y=NULL)
#'
#' print(model)
#' }
#' @section Method Arguments:
#' \describe{
#'   \item{element}{a string with the name of one of the matrices / vectors of the problem}
#'   \item{newPars}{list with optimization parameters}
#' }
#' @details
#' Allows one to solve a parametric
#' problem with for example warm starts between updates of the parameter, c.f. the examples.
#' The object returned by \code{osqp} contains several methods which can be used to either update/get details of the
#' problem, modify the optimization settings or attempt to solve the problem.
#' @return
#' An R6-object of class "osqp_model" with methods defined which can be further
#' used to solve the problem with updated settings / parameters.
#' @examples
#' ## example, adapted from OSQP documentation
#' library(Matrix)
#'
#' P <- Matrix(c(11., 0.,
#'               0., 0.), 2, 2, sparse = TRUE)
#' q <- c(3., 4.)
#' A <- Matrix(c(-1., 0., -1., 2., 3.,
#'               0., -1., -3., 5., 4.)
#'               , 5, 2, sparse = TRUE)
#' u <- c(0., 0., -15., 100., 80)
#' l <- rep_len(-Inf, 5)
#'
#' settings <- osqpSettings(verbose = FALSE)
#'
#' model <- osqp(P, q, A, l, u, settings)
#'
#' # Solve
#' res <- model$Solve()
#'
#' # Define new vector
#' q_new <- c(10., 20.)
#'
#' # Update model and solve again
#' model$Update(q = q_new)
#' res <- model$Solve()
#' @importFrom Matrix sparseMatrix
#' @export
osqp = function(P=NULL, q=NULL, A=NULL, l=NULL, u=NULL, pars = osqpSettings()) {

  if(is.null(P) && is.null(q))
    stop("At least one of P and q must be supplied")

  if (is.null(P))
    n = length(q)
  else {
    dimP <- dim(P)
    n = dim(P)[1L]
    if (dimP[2L] != n) stop("P must be symmetric!")
  }

  if (is.null(P)){
    P = Matrix::sparseMatrix(integer(), integer(), x = numeric(), dims = c(n, n))
  } else {
    ## P = triu(as(P, "dgCMatrix"))
    P <- ensure_dtc_matrix(P)
  }

  if (is.null(q))
    q = numeric(n)
  else
    q = as.numeric(q)


  if (is.null(A)) {
    m = 0
    A = Matrix::sparseMatrix(integer(), integer(), x = numeric(), dims = c(m, n))
    u = l = numeric()
  } else {
    m = nrow(A)
    ## A = as(A, "dgCMatrix")
    A <- ensure_dgc_matrix(A)
    if (is.null(u))
      u = rep_len(Inf, m)
    else
      u = as.numeric(u)

    if (is.null(l))
      l = rep_len(-Inf, m)
    else
      l = as.numeric(l)
  }

  stopifnot(length(q) == n,
            dim(A) == c(m, n),
            length(l) == m,
            length(u) == m,
            l <= u)

  R6Class("osqp_model",
          public =
            list(
              initialize = function(P=NULL, q=NULL, A=NULL, l=NULL, u=NULL, pars=list()) {
                private$.work = osqpSetup(P, q, A, l, u, pars)
              },
              Solve = function() osqpSolve(private$.work),
              Update = function(q = NULL, l = NULL, u = NULL,
                                Px = NULL, Px_idx = NULL, Ax = NULL, Ax_idx = NULL) {
                dims = osqpGetDims(private$.work)
                stopifnot(length(q) %in% c(0, dims[[1]]),
                          length(l) %in% c(0, dims[[2]]),
                          length(u) %in% c(0, dims[[2]])
                          )
                osqpUpdate(private$.work, q, l, u, Px, Px_idx, Ax, Ax_idx)
              },
              GetParams = function() osqpGetParams(private$.work),
              GetDims = function() osqpGetDims(private$.work),
              UpdateSettings = function(newpars = osqpSettings()) {
                stopifnot(is.list(newpars))
                for (i in seq_along(newpars))
                  osqpUpdateSettings(private$.work, newpars[[i]], names(newpars)[[i]])
              },
              GetData = function(element = c("P", "q", "A", "l", "u")) {
                element = match.arg(element)

                osqpGetData(private$.work, element)
              },
              WarmStart = function(x=NULL, y=NULL) {

                dims = osqpGetDims(private$.work)
                stopifnot(length(x) %in% c(0, dims[[1]]),
                          length(y) %in% c(0, dims[[2]]))
                osqpWarmStart(private$.work, x, y)
              }
            ),
          private = list(.work=NULL)
  )$new(P, q, A, l, u, pars)
}

#' @export
format.osqp_model = function(x, ...) {
  dims = x$GetDims()
  sprintf("OSQP-modelobject\n\nNumber of variables: %i\nNumber of constraints: %i", dims[[1]], dims[[2]])
}

#' @export
print.osqp_model = function(x, ...)
  cat(format(x))


private = NULL # to suppress cran note
