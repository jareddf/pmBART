small_fit <- function(x.test = NULL, ...) {
  set.seed(123)
  x <- cbind(x1 = seq(0.02, 0.98, length.out = 20), x2 = runif(20))
  y <- rep(c(0, 1), 10)
  probit_monbart(
    x, y, x.test = x.test,
    ntree = 3, ndpost = 4, nskip = 1, mgsize = 10,
    nkeeptreedraws = 0, ...
  )
}

test_that("a fit has documented probability outputs", {
  test_x <- cbind(x1 = c(0.2, 0.8), x2 = c(0.3, 0.7))
  fit <- small_fit(test_x)

  expect_s3_class(fit, "pmBART_fit")
  expect_identical(dim(fit$prob.train), c(4L, 20L))
  expect_identical(dim(fit$prob.test), c(4L, 2L))
  expect_length(fit$prob.train.mean, 20L)
  expect_length(fit$prob.test.mean, 2L)
  expect_null(fit$yhat.train.mean)
  expect_null(fit$yhat.test.mean)
  expect_true(all(fit$prob.train >= 0 & fit$prob.train <= 1))
  expect_true(all(fit$prob.test >= 0 & fit$prob.test <= 1))
  expect_equal(unname(fit$prob.train.mean), colMeans(fit$prob.train))
  expect_equal(unname(fit$prob.test.mean), colMeans(fit$prob.test))
})

test_that("fitting without test data returns empty test outputs", {
  fit <- small_fit()

  expect_identical(dim(fit$prob.test), c(4L, 0L))
  expect_length(fit$prob.test.mean, 0L)
})

test_that("probability means require retained draws", {
  test_x <- cbind(x1 = c(0.2, 0.8), x2 = c(0.3, 0.7))
  fit <- small_fit(test_x, nkeeptrain = 0, nkeeptest = 0)

  expect_true(all(is.na(fit$prob.train.mean)))
  expect_true(all(is.na(fit$prob.test.mean)))
})

test_that("invalid inputs fail before entering the sampler", {
  x <- matrix(seq_len(20), ncol = 2)

  expect_error(probit_monbart(x, rep(1, 10)), "each class")
  expect_error(probit_monbart(x, c(rep(0, 9), 2)), "only 0 and 1")
  expect_error(
    probit_monbart(x, rep(c(0, 1), 5), x.test = matrix(1, 2, 3)),
    "same number of columns"
  )
  expect_error(
    probit_monbart(x, rep(c(0, 1), 5), ndpost = 5, nkeeptrain = 3),
    "divisor"
  )
})
