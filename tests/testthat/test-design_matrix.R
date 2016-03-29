context("Performing calculations on Design Matrices")

test_that("Polynomial design matrix works fine", {
  polyn <- polynomial.object(M = 2)
  obs <- c(1, 2, 3)
  expect_is(design_matrix(polyn, obs), "list")
  expect_error(design_matrix(polyn, cbind(1, obs)))
  expect_identical(design_matrix(polyn, obs)$H,
                   matrix(c(1, 1, 1, 1, 2, 3, 1, 4, 9), ncol = 3))

  wr_basis <- structure(list(M = 2), class = "basis_x")
  expect_error(design_matrix(wr_basis))
})

test_that("RBF design matrix works fine", {
  rbf <- rbf.object(M = 2, gamma = 1)
  obs <- c(0.1, 0.3, 0.5)
  expect_is(design_matrix(rbf, obs), "list")
  expect_error(design_matrix(rbf, cbind(1, obs)))
  expect_gt(design_matrix(rbf, obs)$basis$mus[1], -0.34)

  wr_basis <- structure(list(M = 2), class = "basis_x")
  expect_error(design_matrix(wr_basis))
})
