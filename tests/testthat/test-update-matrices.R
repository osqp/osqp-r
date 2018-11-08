context("test-update-matrices.R")

library(Matrix)

define_simple_qp <- function(){
    P <- Matrix(c(11., 0.,
                  0., 2.), 2, 2, sparse = TRUE)
    q <- c(3., 4.)
    A <- Matrix(c(-1., 0., -1., 2., 3.,
                  0., -1., -3., 5., 4.)
                  , 5, 2, sparse = TRUE)
    u <- c(0., 0., -15., 100., 80)
    l <- rep_len(-Inf, 5)

    settings <- osqpSettings(verbose = FALSE,
                             eps_rel = 1e-09,
                             eps_abs = 1e-09)

    # Create OSQP model
    model <- osqp(P, q, A, l, u, settings)

    return(model)
}


test_that("Solve update matrices", {

    # Create OSQP model
    model <- define_simple_qp()

    # Solve
    res <- model$Solve()

    expect_equal(res$x, c(0., 5.), 1e-03)
    expect_equal(res$y, c(1.66666, 0., 1.3333, 0., 0.), 1e-03)
    expect_equal(res$info$obj_val, 20., 1e-03)

    # Update P
    Px <- c(1.,  3.)

    model$Update(Px = Px)
    res_P <- model$Solve()

    expect_equal(res$x, c(0., 5.), 1e-03)
    expect_equal(res$y, c(1.66666, 0., 1.3333, 0., 0.), 1e-03)
    expect_equal(res$info$obj_val, 20., 1e-03)

})

