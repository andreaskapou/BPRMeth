context("Performing calculations of BPR likelihood function")

test_that("bpr_likelihood function works fine", {
  polyn <- create_polynomial_object(M = 1)
  obs <- c(0, .5, .7)
  des_mat <- .design_matrix(polyn, obs)
  data <- matrix(c(10, 12, 15, 7, 9, 8), ncol = 2)
  w <- c(.1, .1)
  expect_gt(.bpr_likelihood(w, des_mat$H, data), -5.8)
  expect_lt(.bpr_likelihood(w, des_mat$H, data), -5.7)
})

test_that("bpr_gradient function works fine", {
  obj <- create_rbf_object(M = 2)
  obs <- c(0, .2, .5)
  des_mat <- .design_matrix(obj, obs)
  data <- matrix(c(10, 12, 15, 7, 9, 8), ncol = 2)
  w <- c(.1, .1, 1)
  expect_gt(.bpr_gradient(w, des_mat$H, data)[1], -15.4)
  expect_lt(.bpr_gradient(w, des_mat$H, data)[1], -15.3)
})
