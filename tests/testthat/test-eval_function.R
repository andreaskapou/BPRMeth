context("Checking for correct evaluation of basis functions")

test_that("polynomial evaluation function works fine", {
  x <- create_polynomial_object(M=2)
  obs <- c(1,2,3)
  w <- c(1, 1, 1)
  expect_identical(eval_function(x, obs, w), c(3, 7, 13))
  expect_error(eval_function.rbf(x, obs, w))
})

test_that("rbf evaluation function works fine", {
  x <- create_rbf_object(M=2)
  obs <- c(1,2,3)
  w <- c(1, 1, 1)
  x$mus <- c(0, 0)
  expect_gt(eval_function(x, obs, w)[1], 1.73)
  expect_lt(eval_function(x, obs, w)[1], 1.74)
  expect_error(eval_function.polynomial(x, obs, w))
})
