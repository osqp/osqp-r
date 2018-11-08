context("test-solve-basic-qp.R")

library(Matrix)

define_simple_qp <- function(){
    P <- Matrix(c(11., 0.,
                  0., 0.), 2, 2, sparse = TRUE)
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


test_that("Solve basic QP", {

    # Create OSQP model
    model <- define_simple_qp()

    # Solve
    res <- model$Solve()

    expect_equal(res$x, c(0., 5.), 1e-03)
    expect_equal(res$y, c(1.66666, 0., 1.3333, 0., 0.), 1e-03)
    expect_equal(res$info$obj_val, 20., 1e-03)

})

test_that("Update basic QP", {

    # Create OSQP model
    model <- define_simple_qp()

    # Solve
    res <- model$Solve()

    # Define new vector
    q_new <- c(10., 20.)

    # Update model and solve again
    model$Update(q = q_new)
    res <- model$Solve()

    expect_equal(res$x, c(0., 5.), 1e-03)
    expect_equal(res$y, c(3.333, 0., 6.666, 0., 0.), 1e-03)
    expect_equal(res$info$obj_val, 100., 1e-03)

})

test_that("Update bounds QP", {

    # Create OSQP model
    model <- define_simple_qp()

    # Solve
    res <- model$Solve()

    # Define new vector
    l_new <- -100 * rep(1, 5)
    u_new <- 1000 * rep(1, 5)

    # Update model and solve again
    model$Update(l = l_new, u = u_new)
    res_updated <- model$Solve()

    expect_equal(res_updated$x, c(-0.12727273,
                                  -19.94909091), 1e-03)
    expect_equal(res_updated$y, c(0., 0., 0.,
                                  -0.8, 0.), 1e-03)
    expect_equal(res_updated$info$obj_val,
                 -80.0890909023583, 1e-03)

})

