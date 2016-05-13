context("Performing calculations on Basis functions")

test_that("polynomial basis object works fine", {
  expect_error(create_polynomial_object(1.5))
  expect_error(create_polynomial_object(-1))
  expect_error(create_polynomial_object("non-numeric"))
  expect_is(create_polynomial_object(2), "polynomial")
})

test_that("radial basis object works fine", {
  expect_error(create_rbf_object(1.5))
  expect_error(create_rbf_object(-1))
  expect_error(create_rbf_object("non-numeric"))
  expect_error(create_rbf_object(eq_spaced_mus = "non-logical"))
  expect_error(create_rbf_object(gamma = "non-numeric"))
  expect_error(create_rbf_object(M = 2, gamma = 0))
  expect_error(create_rbf_object(M = 2, gamma = 0, mus = c(2,3,4)))
  expect_is(create_rbf_object(M = 2, mus = c(2,2)), "rbf")
  expect_is(create_rbf_object(2), "rbf")
})

test_that("polynomial function works fine", {
  expect_identical(.polynomial_basis(c(4), 2), 16)
  expect_identical(.polynomial_basis(c(2,6), 2), c(4, 36))
  expect_identical(.polynomial_basis(matrix(c(1,2,3,4), ncol = 2), 2),
                   matrix(c(1, 4, 9, 16), ncol = 2))
})

test_that("rbf function works fine", {
  expect_identical(.rbf_basis(1,1,1), 1)
  expect_gt(.rbf_basis(1,2,1), 0.366)
  expect_gt(.rbf_basis(1,5,1), 1.1e-07)
  expect_gt(.rbf_basis(1,5,.1), 0.2)
  expect_gt(.rbf_basis(c(2,2), c(2,1), 1), 0.366)
  expect_gt(.rbf_basis(c(2,3), c(4,4), 1), 0.0067)
})
