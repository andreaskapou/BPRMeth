context("Performing calculations of BPR likelihood function")

test_that("bpr_likelihood function works fine", {
  polyn <- polynomial.object(M = 1)
  obs <- c(0,.5,.7)
  des_mat <- design_matrix(polyn, obs)
  data <- matrix(c(10,12,15,7,9,8), ncol=2)
  w <- c(.1,.1)
  expect_gt(bpr_likelihood(w, des_mat$H, data), -5.8)
  expect_lt(bpr_likelihood(w, des_mat$H, data), -5.7)
})
