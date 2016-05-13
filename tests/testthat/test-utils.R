context("Checking utils functions work")

test_that("data centering works", {
  region <- c(2, 4, 5)
  tss <- 3
  expect_identical(.center_loc(region, tss, "+"), c(-1, 1, 2))
  expect_identical(.center_loc(region, tss, "-"), c(1, -1, -2))
  expect_error(.center_loc(region, tss, 9))
  expect_error(.center_loc(region, "tss", 9))
  expect_identical(.center_loc(21, 2, "+"), 19)
})

test_that("minmax scaling works",{
  data <- c(-20, 0, 15, 20)
  xmin <- -20
  xmax <- 20
  expect_identical(.minmax_scaling(data, xmin, xmax, fmin=0, fmax=1), c(0, 0.5, 0.875, 1))
  data <- c(-20, 20)
  expect_identical(.minmax_scaling(data, fmin=-2, fmax=2), c(-2, 2))
  expect_identical(.minmax_scaling(data, -xmin, xmax, fmin=0, fmax=1), c(-20, 20))
  expect_identical(.minmax_scaling(4, 0, 8, fmin=0, fmax=1), 0.5)
})

test_that("partition data works",{
  set.seed(123)
  x <- matrix(c(-20, 0, 15, 20, 20, 10, 4, 5), ncol=2)
  y <- c(2, 5, 0, 8)
  expect_identical(.partition_data(x, y, c(4,1,3))$train$y, c(8,2,0))
  expect_identical(.partition_data(x, y, c(4,1,3))$test$y, c(5))
  expect_identical(.partition_data(x, y, c(4,1,3))$test$x.2, c(10))
})

test_that("calculate errors works",{
  x <- c(4, 6, 9, 10, 4, 6, 4, 7, 8, 7)
  y <- c(5, 6, 8, 10, 4, 8, 4, 9, 8, 9)
  expect_identical(.calculate_errors(x, y)$mae, 0.8)
  expect_identical(.calculate_errors(x, y)$mse, 1.4)
  expect_gt(.calculate_errors(x, y)$pcc, 0.87)
  expect_lt(.calculate_errors(x, y)$pcc, 0.88)
})


test_that("log sum exp function works",{
  x <- c(0.001, 0.5, 2, 1.4, 1.5)
  expect_gt(.log_sum_exp(x), 2.92)
  expect_lt(.log_sum_exp(x), 2.93)
})
