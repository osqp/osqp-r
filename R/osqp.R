
#' OSQP Model Class (S7)
#'
#' @description
#' An S7 class representing an OSQP solver model. Methods are accessed via
#' computed properties using the `@` operator.
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom Matrix triu
#' @importFrom methods as
#' @import S7
#' @keywords internal
OSQP_Model <- new_class("OSQP_Model",
  package = NULL,
  properties = list(
    .work = new_property(class_any),

    Solve = new_property(class_function, getter = function(self) {
      work <- self@.work
      function() osqpSolve(work)
    }),

    Update = new_property(class_function, getter = function(self) {
      work <- self@.work
      function(q = NULL, l = NULL, u = NULL,
               Px = NULL, Px_idx = NULL, Ax = NULL, Ax_idx = NULL) {
        dims <- osqpGetDims(work)
        n <- dims[[1]]; m <- dims[[2]]
        if (!length(q) %in% c(0L, n))
          stop(cli::format_error("{.arg q} must have length {n} (number of variables), not {length(q)}."))
        if (!length(l) %in% c(0L, m))
          stop(cli::format_error("{.arg l} must have length {m} (number of constraints), not {length(l)}."))
        if (!length(u) %in% c(0L, m))
          stop(cli::format_error("{.arg u} must have length {m} (number of constraints), not {length(u)}."))
        osqpUpdate(work, q, l, u, Px, Px_idx, Ax, Ax_idx)
      }
    }),

    WarmStart = new_property(class_function, getter = function(self) {
      work <- self@.work
      function(x = NULL, y = NULL) {
        dims <- osqpGetDims(work)
        n <- dims[[1]]; m <- dims[[2]]
        if (!length(x) %in% c(0L, n))
          stop(cli::format_error("{.arg x} must have length {n} (number of variables), not {length(x)}."))
        if (!length(y) %in% c(0L, m))
          stop(cli::format_error("{.arg y} must have length {m} (number of constraints), not {length(y)}."))
        osqpWarmStart(work, x, y)
      }
    }),

    ColdStart = new_property(class_function, getter = function(self) {
      work <- self@.work
      function() osqpColdStart(work)
    }),

    GetParams = new_property(class_function, getter = function(self) {
      work <- self@.work
      function() osqpGetParams(work)
    }),

    GetDims = new_property(class_function, getter = function(self) {
      work <- self@.work
      function() osqpGetDims(work)
    }),

    UpdateSettings = new_property(class_function, getter = function(self) {
      work <- self@.work
      function(newpars = osqpSettings()) {
        if (!is.list(newpars))
          stop(cli::format_error("{.arg newpars} must be a list, not {.cls {class(newpars)}}."))
        for (i in seq_along(newpars))
          osqpUpdateSettings(work, newpars[[i]], names(newpars)[[i]])
      }
    }),

    GetData = new_property(class_function, getter = function(self) {
      work <- self@.work
      function(element = c("P", "q", "A", "l", "u")) {
        element <- match.arg(element)
        osqpGetData(work, element)
      }
    })
  ),
  constructor = function(work) {
    new_object(S7_object(), .work = work)
  }
)


#' OSQP Solver object
#'
#' @param P,A sparse matrices of class dgCMatrix or coercible into such, with P positive semidefinite. (In the interest of efficiency, only the upper triangular part of P is used)
#' @param q,l,u Numeric vectors, with possibly infinite elements in l and u
#' @param pars list with optimization parameters, conveniently set with the function
#' \code{\link{osqpSettings}}. For \code{model@@UpdateSettings(newPars)} only a subset of the settings
#' can be updated once the problem has been initialized.
#' @seealso \code{\link{solve_osqp}}
#' @section Usage:
#' \preformatted{model = osqp(P=NULL, q=NULL, A=NULL, l=NULL, u=NULL, pars=osqpSettings())
#'
#' model@@Solve()
#' model@@Update(q = NULL, l = NULL, u = NULL, Px = NULL, Px_idx = NULL, Ax = NULL, Ax_idx = NULL)
#' model@@GetParams()
#' model@@GetDims()
#' model@@UpdateSettings(newPars = list())
#'
#' model@@GetData(element = c("P", "q", "A", "l", "u"))
#' model@@WarmStart(x=NULL, y=NULL)
#' model@@ColdStart()
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
#' The object returned by \code{osqp} contains several computed properties (accessed via \code{@@})
#' which can be used to either update/get details of the
#' problem, modify the optimization settings or attempt to solve the problem.
#' @return
#' An S7 object of class "OSQP_Model" with computed properties that return methods.
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
#' res <- model@@Solve()
#'
#' # Define new vector
#' q_new <- c(10., 20.)
#'
#' # Update model and solve again
#' model@@Update(q = q_new)
#' res <- model@@Solve()
#' @export
osqp <- function(P = NULL, q = NULL, A = NULL, l = NULL, u = NULL, pars = osqpSettings()) {

  if (is.null(P) && is.null(q))
    stop(cli::format_error("At least one of {.arg P} and {.arg q} must be supplied."))

  if (is.null(P))
    n <- length(q)
  else {
    dimP <- dim(P)
    n <- dimP[1L]
    if (dimP[2L] != n) stop(cli::format_error("{.arg P} must be a square matrix, got {dimP[1L]} x {dimP[2L]}."))
  }

  if (is.null(P)) {
    P <- Matrix::sparseMatrix(integer(), integer(), x = numeric(), dims = c(n, n))
  } else {
    P <- ensure_dtc_matrix(P)
  }

  if (is.null(q))
    q <- numeric(n)
  else
    q <- as.numeric(q)

  if (is.null(A)) {
    m <- 0L
    A <- Matrix::sparseMatrix(integer(), integer(), x = numeric(), dims = c(m, n))
    u <- l <- numeric()
  } else {
    m <- nrow(A)
    A <- ensure_dgc_matrix(A)
    if (is.null(u))
      u <- rep_len(Inf, m)
    else
      u <- as.numeric(u)

    if (is.null(l))
      l <- rep_len(-Inf, m)
    else
      l <- as.numeric(l)
  }

  if (length(q) != n)
    stop(cli::format_error("{.arg q} must have length {n}, not {length(q)}."))
  if (!identical(dim(A), c(m, n)))
    stop(cli::format_error("{.arg A} must be {m} x {n}, got {nrow(A)} x {ncol(A)}."))
  if (length(l) != m)
    stop(cli::format_error("{.arg l} must have length {m}, not {length(l)}."))
  if (length(u) != m)
    stop(cli::format_error("{.arg u} must have length {m}, not {length(u)}."))
  if (any(l > u))
    stop(cli::format_error("Lower bounds {.arg l} must not exceed upper bounds {.arg u}."))

  work <- osqpSetup(P, q, A, l, u, pars)
  OSQP_Model(work)
}

#' @export
format.OSQP_Model <- function(x, ...) {
  dims <- x@GetDims()
  sprintf("OSQP-model object\n\nNumber of variables: %i\nNumber of constraints: %i", dims[[1]], dims[[2]])
}

#' @export
print.OSQP_Model <- function(x, ...) {
  dims <- x@GetDims()
  cli::cli_h3("OSQP Model")
  cli::cli_ul(c(
    "Number of variables: {dims[[1]]}",
    "Number of constraints: {dims[[2]]}"
  ))
  invisible(x)
}

## ------------------------------------------------------------------
## Deprecation bridge: allow model$Method() with a warning.
## TODO: Remove this bridge and tests/testthat/test-dollar-deprecation.R
##       after the deprecation period (target: next CRAN release after 1.0.0).
## ------------------------------------------------------------------
.osqp_methods <- c("Solve", "Update", "WarmStart", "ColdStart",
                    "GetParams", "GetDims", "UpdateSettings", "GetData")

#' @export
`$.OSQP_Model` <- function(x, name) {
  if (name %in% .osqp_methods) {
    warning(cli::format_warning(
      c("!" = "Using {.code model${name}()} is deprecated for {.cls OSQP_Model} objects.",
        "i" = "Use {.code model@{name}()} instead.")),
      call. = FALSE
    )
    return(prop(x, name))
  }
  NULL
}
