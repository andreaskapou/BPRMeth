context("Performing calculations of likelihood functions")

test_that("bpr_likelihood function works", {
    polyn <- create_polynomial_object(M = 1)
    obs <- c(0, .5, .7)
    des_mat <- design_matrix(polyn, obs)
    X <- matrix(c(10, 12, 15, 7, 9, 8), ncol = 2)
    X <- cbind(obs, X)
    w <- c(.1, .1)
    expect_gt(bpr_log_likelihood(w, X, des_mat$H, 1/2, FALSE), -5.8)
    expect_lt(bpr_log_likelihood(w, X, des_mat$H, 1/2, FALSE), -5.7)
})

test_that("bpr_gradient function works", {
    obj <- create_rbf_object(M = 2)
    obs <- c(0, .2, .5)
    des_mat <- design_matrix(obj, obs)
    X <- matrix(c(10, 12, 15, 7, 9, 8), ncol = 2)
    X <- cbind(obs, X)
    w <- c(.1, .1, 1)
    expect_gt(bpr_gradient(w, X, des_mat$H, 1/2, FALSE)[1], -15.4)
    expect_lt(bpr_gradient(w, X, des_mat$H, 1/2, FALSE)[1], -15.2)
})


test_that("betareg_likelihood function works", {
    obj <- create_rbf_object(M = 2)
    obs <- c(0, .2, .5)
    H <- design_matrix(obj, obs)$H
    X <- matrix(c(0.1, 0.5, 0.3, 4, 4, 4), ncol = 2)
    X <- cbind(obs, X)
    w <- c(.1, .6, 2)
    expect_gt(betareg_log_likelihood(w, X, H, 1/2, FALSE), -23.9)
    expect_lt(betareg_log_likelihood(w, X, H, 1/2, FALSE), -23.8)
})


test_that("betareg_gradient function works", {
    obj <- create_rbf_object(M = 2)
    obs <- c(0, .2, .5)
    H <- design_matrix(obj, obs)$H
    X <- matrix(c(0.1, 0.5, 0.3, 4, 4, 4), ncol = 2)
    X <- cbind(obs, X)
    w <- c(.1, .6, 2)
    expect_gt(betareg_gradient(w, X, H, 1/2, FALSE)[1], -9.1)
    expect_lt(betareg_gradient(w, X, H, 1/2, FALSE)[1], -9.05)
})
