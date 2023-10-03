

context("Parameter settings")

test_that("Parameter settings and types work", {
  ## https://github.com/osqp/osqp-r/issues/17
  max_iter  <-  1000.0
  eps  <- 1e-5
  settings  <- osqp::osqpSettings(max_iter = max_iter, eps_abs = eps)
  expect_identical(settings, list(max_iter = 1000L, eps_abs = 1e-5))
})


test_that("Lower bounds can't exceed upper bounds etc.", {
  ## https://github.com/osqp/osqp-r/issues/29
  P <- Matrix::Matrix(c(11., 0., 0., 0.), 2, 2, sparse = TRUE)
  q <- c(3., 4.)
  A <- Matrix::Matrix(c(-1., 0., -1., 2., 3., 0., -1., -3., 5., 4.) , 5, 2, sparse = TRUE)
  u <- c(0., 0., -15., 100., 80)
  l <- u + 1
  settings <- osqpSettings(verbose = TRUE)
  # Solve with OSQP.
  ## Should throw error, but not crash as it is now caught at the R level, not deep in C++
  expect_error(solve_osqp(P, q, A, l, u, settings))
})


test_that("Time limit setting works", {
  P <- 2 * matrix(c(0.25^2, .04, -.005, .04, .45^2, -.01, -.005, -.01, .05^2), nrow = 3)
  A <- matrix(c(1, 21, 1, 0, 0, 1, 30, 0, 1, 0, 1, 8, 0, 0, 1), nrow = 5)
  l <- c(1, 18, 0, 0, 0)
  u <- c(1, rep(Inf, 4))
  res <- solve_osqp(P = P, A = A, l = l, u = u, pars = osqpSettings(time_limit = 1e-10))
  expect_equal(res$info$status, "run time limit reached")  ## runtime limit reached
})


