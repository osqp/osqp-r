
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
  P <- list(i = c(1L, 2L, 1L, 2L), j = c(1L, 1L, 2L, 2L), v = c(11.0, 2.0, 2.0, 1.0))
  P$nrow <- 2L; P$ncol <- 2L; class(P) <- "simple_triplet_matrix"
  q <- c(3., 4.)
  A <- list(i = c(1L, 3L, 4L, 5L, 2L, 3L, 4L, 5L),
            j = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L),
            v = c(-1.0, -1.0, 2.0, 3.0, -1.0, -3.0, 5.0, 4.0))
  A$nrow <- 5L; A$ncol <- 2L; class(A) <- "simple_triplet_matrix"
  u <- c(0., 0., -15., 100., 80)
  l <- rep_len(-Inf, 5)

  settings <- osqpSettings(verbose = TRUE)
  res <- solve_osqp(P, q, A, l, u, settings)
  expect_equal(res$info$status, "solved")
})

