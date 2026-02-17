## ---------------------------------------------------------------
## DEPRECATION BRIDGE TESTS: $ (old R6 API) -> @ (S7 API)
##
## These tests ensure that code written for the old R6-based osqp
## (which used model$Method()) continues to work via the $.OSQP_Model
## bridge defined in R/osqp.R.
##
## TODO: Remove this file and $.OSQP_Model in R/osqp.R after the
##       deprecation period (target: next CRAN release after 1.0.0).
## ---------------------------------------------------------------

context("test-dollar-deprecation.R")

library(Matrix)

make_model <- function() {
  P <- Matrix(c(11., 0.,
                0., 0.), 2, 2, sparse = TRUE)
  q <- c(3., 4.)
  A <- Matrix(c(-1., 0., -1., 2., 3.,
                0., -1., -3., 5., 4.),
              5, 2, sparse = TRUE)
  u <- c(0., 0., -15., 100., 80)
  l <- rep_len(-Inf, 5)
  osqp(P, q, A, l, u, osqpSettings(verbose = FALSE))
}

## -- Basic dispatch and warning -----------------------------------

test_that("$ access works but emits deprecation warning", {
  model <- make_model()

  res_at <- model@Solve()
  expect_warning(res_dollar <- model$Solve(), "deprecated")
  expect_equal(res_dollar$x, res_at$x, tolerance = 1e-3)
})

test_that("$ deprecation warning for each method", {
  model <- make_model()

  expect_warning(model$GetDims(), "deprecated")
  expect_warning(model$GetParams(), "deprecated")
  expect_warning(model$GetData("q"), "deprecated")
  expect_warning(model$Update(q = c(3., 4.)), "deprecated")
  expect_warning(model$WarmStart(x = c(0, 0)), "deprecated")
  expect_warning(model$ColdStart(), "deprecated")
  expect_warning(model$UpdateSettings(list(verbose = FALSE)), "deprecated")
})

test_that("$ on unknown name returns NULL", {
  model <- make_model()
  expect_null(model$nonexistent)
})

## -- Master-branch example: solve, update, re-solve ---------------

test_that("master osqp.R example works via $", {
  P <- Matrix(c(11., 0.,
                0., 0.), 2, 2, sparse = TRUE)
  q <- c(3., 4.)
  A <- Matrix(c(-1., 0., -1., 2., 3.,
                0., -1., -3., 5., 4.),
              5, 2, sparse = TRUE)
  u <- c(0., 0., -15., 100., 80)
  l <- rep_len(-Inf, 5)

  model <- osqp(P, q, A, l, u, osqpSettings(verbose = FALSE))

  suppressWarnings(res <- model$Solve())
  expect_equal(res$info$status, "solved")

  q_new <- c(10., 20.)
  suppressWarnings(model$Update(q = q_new))
  suppressWarnings(res2 <- model$Solve())
  expect_equal(res2$info$status, "solved")
})

## -- GetData for every element ------------------------------------

test_that("$ GetData returns all elements", {
  model <- make_model()
  for (el in c("P", "q", "A", "l", "u")) {
    suppressWarnings(val <- model$GetData(el))
    expect_true(length(val) > 0, info = paste("GetData element:", el))
  }
})

## -- UpdateSettings round-trip ------------------------------------

test_that("$ UpdateSettings + GetParams round-trip", {
  model <- make_model()
  suppressWarnings(model$UpdateSettings(osqpSettings(max_iter = 2000L)))
  suppressWarnings(pars <- model$GetParams())
  expect_equal(pars$max_iter, 2000L)
})

## -- WarmStart / ColdStart cycle ----------------------------------

test_that("$ WarmStart and ColdStart cycle", {
  model <- make_model()
  suppressWarnings(res <- model$Solve())
  suppressWarnings(model$WarmStart(x = res$x, y = res$y))
  suppressWarnings(res2 <- model$Solve())
  expect_equal(res2$info$status, "solved")

  suppressWarnings(model$ColdStart())
  suppressWarnings(res3 <- model$Solve())
  expect_equal(res3$info$status, "solved")
})

## -- Matrix update (Px, Ax) via $ ---------------------------------

test_that("$ Update with Px and Ax", {
  P <- Matrix(c(4., 1.,
                1., 2.), 2, 2, sparse = TRUE)
  q <- c(1., 1.)
  A <- Matrix(c(1., 1., 0.,
                1., 0., 1.), 3, 2, sparse = TRUE)
  l <- c(1., 0., 0.)
  u <- c(1., 0.7, 0.7)

  model <- osqp(P, q, A, l, u, osqpSettings(verbose = FALSE))
  suppressWarnings(res1 <- model$Solve())
  expect_equal(res1$info$status, "solved")

  P_new <- Matrix(c(5., 1.5,
                    1.5, 1.), 2, 2, sparse = TRUE)
  A_new <- Matrix(c(1.2, 1.5, 0.,
                    1.1, 0., 0.8), 3, 2, sparse = TRUE)

  suppressWarnings(model$Update(Px = triu(P_new)@x, Ax = A_new@x))
  suppressWarnings(res2 <- model$Solve())
  expect_equal(res2$info$status, "solved")
})
