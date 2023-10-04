
context("Sparse Coercions")

test_that("Dense P & A are properly coerced", {
  P <- matrix(c(11., 2.0, 2.0, 1.0), 2, 2)
  q <- c(3., 4.)
  A <- matrix(c(-1., 0., -1., 2., 3., 0., -1., -3., 5., 4.) , 5, 2)
  u <- c(0., 0., -15., 100., 80)
  l <- rep_len(-Inf, 5)

  settings <- osqpSettings(verbose = TRUE)
  res <- solve_osqp(P, q, A, l, u, settings)
  expect_equal(res$info$status, "solved")
})

test_that("Triplet-form P & A are properly coerced", {
  P <- slam::as.simple_triplet_matrix(matrix(c(11., 2.0, 2.0, 1.0), 2, 2))
  q <- c(3., 4.)
  A <- slam::as.simple_triplet_matrix(matrix(c(-1., 0., -1., 2., 3., 0., -1., -3., 5., 4.) , 5, 2))
  u <- c(0., 0., -15., 100., 80)
  l <- rep_len(-Inf, 5)

  settings <- osqpSettings(verbose = TRUE)
  res <- solve_osqp(P, q, A, l, u, settings)
  expect_equal(res$info$status, "solved")
})

