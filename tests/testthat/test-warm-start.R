context("Warm start and cold start")

library(Matrix)

define_simple_qp <- function() {
  P <- Matrix(c(11., 0.,
                0., 0.), 2, 2, sparse = TRUE)
  q <- c(3., 4.)
  A <- Matrix(c(-1., 0., -1., 2., 3.,
                0., -1., -3., 5., 4.),
              5, 2, sparse = TRUE)
  u <- c(0., 0., -15., 100., 80)
  l <- rep_len(-Inf, 5)

  settings <- osqpSettings(verbose = FALSE,
                           eps_rel = 1e-05,
                           eps_abs = 1e-05,
                           warm_starting = TRUE)

  model <- osqp(P, q, A, l, u, settings)
  model
}

test_that("WarmStart with primal variables speeds up re-solve", {
  model <- define_simple_qp()

  # Cold solve
  res1 <- model@Solve()
  expect_equal(res1$info$status, "solved")
  iter_cold <- res1$info$iter

  # Warm start with the solution and re-solve
  model@WarmStart(x = res1$x, y = res1$y)
  res2 <- model@Solve()
  expect_equal(res2$info$status, "solved")
  expect_equal(res2$x, res1$x, tolerance = 1e-03)

  # Warm-started solve should take fewer iterations
  expect_lte(res2$info$iter, iter_cold)
})

test_that("WarmStart with only x works", {
  model <- define_simple_qp()
  res1 <- model@Solve()

  model@WarmStart(x = res1$x)
  res2 <- model@Solve()
  expect_equal(res2$info$status, "solved")
  expect_equal(res2$x, res1$x, tolerance = 1e-03)
})

test_that("WarmStart with only y works", {
  model <- define_simple_qp()
  res1 <- model@Solve()

  model@WarmStart(y = res1$y)
  res2 <- model@Solve()
  expect_equal(res2$info$status, "solved")
  expect_equal(res2$x, res1$x, tolerance = 1e-03)
})

test_that("ColdStart resets iterate and re-solves correctly", {
  model <- define_simple_qp()

  res1 <- model@Solve()
  expect_equal(res1$info$status, "solved")

  # Warm start, then cold start to reset
  model@WarmStart(x = res1$x, y = res1$y)
  model@ColdStart()
  res2 <- model@Solve()

  expect_equal(res2$info$status, "solved")
  expect_equal(res2$x, res1$x, tolerance = 1e-03)
})

test_that("WarmStart after parameter update works", {
  model <- define_simple_qp()

  # Solve original problem
  res1 <- model@Solve()

  # Update q and warm start with previous solution
  q_new <- c(10., 20.)
  model@Update(q = q_new)
  model@WarmStart(x = res1$x, y = res1$y)

  res2 <- model@Solve()
  expect_equal(res2$info$status, "solved")
  expect_equal(res2$x, c(0., 5.), tolerance = 1e-03)
})

test_that("WarmStart rejects wrong-length inputs", {
  model <- define_simple_qp()

  expect_error(model@WarmStart(x = c(1.0)),
               "must have length 2")
  expect_error(model@WarmStart(y = c(1.0, 2.0)),
               "must have length 5")
})
