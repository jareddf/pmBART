#' Fit a probit monotone BART model
#'
#' Fits a Bayesian additive regression tree model for a binary outcome. Every
#' predictor is constrained to have a non-decreasing effect on the latent
#' probit-scale response.
#'
#' @param x.train Numeric matrix or data frame of training predictors. Rows are
#'   observations and columns are predictors.
#' @param y.train Binary response vector. Values must be `0` or `1`.
#' @param x.test Optional numeric matrix or data frame of predictors at which to
#'   save posterior predictions. It must have the same columns, in the same
#'   order, as `x.train`.
#' @param sigest,sigdf,sigquant,sigmaf,lambda Legacy compatibility arguments.
#'   They are currently ignored because the latent probit variance is fixed at
#'   one.
#' @param k Prior shrinkage parameter. Larger values shrink individual trees
#'   more strongly.
#' @param power,base Parameters controlling the tree-depth prior.
#' @param offset Optional probit-scale intercept. By default it is
#'   `qnorm(mean(y.train))`.
#' @param ntree Number of trees in the ensemble.
#' @param ndpost Number of posterior draws after burn-in.
#' @param nskip Number of burn-in iterations.
#' @param mgsize Number of support points in the discrete leaf-parameter prior.
#' @param nkeeptrain Number of posterior training draws to retain. Use zero to
#'   retain none. The value must divide `ndpost`.
#' @param nkeeptest Number of posterior test draws to retain. Use zero to retain
#'   none. The value must divide `ndpost`.
#' @param nkeeptestmean Number of test iterations used to calculate the native
#'   posterior mean. Use zero to disable it. The value must divide `ndpost`.
#' @param nkeeptreedraws Number of serialized tree draws to retain. Use zero to
#'   retain none. The value must divide `ndpost`.
#' @param printevery Print sampler progress every this many iterations. Use zero
#'   (the default) for quiet operation.
#'
#' @return An object of class `"pmBART_fit"` containing posterior latent-scale
#'   fits and probabilities. The most commonly used components are:
#'   \describe{
#'     \item{prob.train}{Matrix of retained training probability draws.}
#'     \item{prob.train.mean}{Posterior mean training probabilities.}
#'     \item{prob.test}{Matrix of retained test probability draws.}
#'     \item{prob.test.mean}{Posterior mean test probabilities.}
#'     \item{yhat.train, yhat.test}{Corresponding latent probit-scale draws.}
#'     \item{treedraws}{Serialized trees and predictor cutpoints.}
#'   }
#'
#' @details
#' The monotonicity constraint is increasing for every predictor. Reverse a
#' predictor before fitting if scientific knowledge requires its effect to be
#' decreasing. Missing or infinite values are not supported.
#'
#' `nkeeptrain`, `nkeeptest`, `nkeeptestmean`, and `nkeeptreedraws` control
#' thinning independently. Each nonzero value must be no larger than, and
#' divide, `ndpost`.
#'
#' @examples
#' set.seed(42)
#' x <- cbind(exposure = runif(40), score = runif(40))
#' p <- pnorm(-1 + 1.5 * x[, "exposure"] + x[, "score"])
#' y <- rbinom(40, size = 1, prob = p)
#'
#' fit <- probit_monbart(
#'   x[1:30, ], y[1:30],
#'   x.test = x[31:40, ],
#'   ntree = 10, ndpost = 20, nskip = 5,
#'   nkeeptreedraws = 0
#' )
#' fit
#' head(fit$prob.test.mean)
#'
#' @export
probit_monbart <- function(
    x.train, y.train, x.test = NULL,
    sigest = NA, sigdf = 3, sigquant = 0.90,
    k = 2,
    power = 0.8, base = 0.25,
    sigmaf = NA,
    lambda = NA,
    offset = NULL,
    ntree = 200,
    ndpost = 1000, nskip = 100,
    mgsize = 50,
    nkeeptrain = ndpost, nkeeptest = ndpost,
    nkeeptestmean = ndpost,
    nkeeptreedraws = ndpost,
    printevery = 0
) {
  x.train <- .validate_predictors(x.train, "x.train")
  n <- nrow(x.train)
  p <- ncol(x.train)

  if (length(y.train) != n) {
    stop("`y.train` must have one value for each row of `x.train`.", call. = FALSE)
  }
  if (anyNA(y.train) || !all(y.train %in% c(0, 1))) {
    stop("`y.train` must contain only 0 and 1, with no missing values.", call. = FALSE)
  }
  if (length(unique(y.train)) < 2L) {
    stop("`y.train` must contain at least one observation from each class.", call. = FALSE)
  }
  y.train <- as.numeric(y.train)

  if (is.null(x.test)) {
    x.test <- matrix(numeric(), nrow = 0L, ncol = p)
    colnames(x.test) <- colnames(x.train)
  } else {
    x.test <- .validate_predictors(x.test, "x.test", allow_zero_rows = TRUE)
    if (ncol(x.test) != p) {
      stop("`x.test` must have the same number of columns as `x.train`.", call. = FALSE)
    }
    train_names <- colnames(x.train)
    test_names <- colnames(x.test)
    if (!is.null(train_names) && !is.null(test_names) &&
        !identical(train_names, test_names)) {
      stop("The column names and order of `x.test` must match `x.train`.", call. = FALSE)
    }
  }

  .check_positive_integer(ntree, "ntree")
  .check_positive_integer(ndpost, "ndpost")
  .check_nonnegative_integer(nskip, "nskip")
  .check_positive_integer(mgsize, "mgsize")
  .check_nonnegative_integer(printevery, "printevery")

  keep <- c(
    nkeeptrain = nkeeptrain,
    nkeeptest = nkeeptest,
    nkeeptestmean = nkeeptestmean,
    nkeeptreedraws = nkeeptreedraws
  )
  for (nm in names(keep)) {
    .check_keep(keep[[nm]], nm, ndpost)
  }

  .check_scalar(k, "k", lower = 0, open_lower = TRUE)
  .check_scalar(base, "base", lower = 0, upper = 1, open_lower = TRUE)
  .check_scalar(power, "power", lower = 0)
  if (is.null(offset)) {
    offset <- stats::qnorm(mean(y.train))
  }
  .check_scalar(offset, "offset")

  tau <- sqrt(1.467) * 3 / (k * sqrt(ntree))
  nu <- 3
  lambda <- 1

  res <- cmonbart(
    t(x.train), y.train, t(x.test),
    tau, nu, lambda, base, power, offset,
    ndpost, nskip, ntree, mgsize,
    nkeeptrain, nkeeptest, nkeeptestmean,
    nkeeptreedraws, printevery
  )

  res$yhat.train <- res$yhat.train + offset
  res$yhat.train.mean <- res$yhat.train.mean + offset
  res$prob.train <- stats::pnorm(res$yhat.train)
  res$prob.train.mean <- if (nrow(res$prob.train) > 0L) {
    colMeans(res$prob.train)
  } else {
    stats::pnorm(res$yhat.train.mean)
  }

  res$yhat.test <- res$yhat.test + offset
  res$yhat.test.mean <- res$yhat.test.mean + offset
  res$prob.test <- stats::pnorm(res$yhat.test)
  res$prob.test.mean <- if (nrow(res$prob.test) > 0L) {
    colMeans(res$prob.test)
  } else if (nkeeptestmean > 0L) {
    stats::pnorm(res$yhat.test.mean)
  } else {
    rep(NA_real_, nrow(x.test))
  }

  colnames(res$yhat.train) <- rownames(x.train)
  colnames(res$prob.train) <- rownames(x.train)
  names(res$yhat.train.mean) <- rownames(x.train)
  names(res$prob.train.mean) <- rownames(x.train)
  colnames(res$yhat.test) <- rownames(x.test)
  colnames(res$prob.test) <- rownames(x.test)
  names(res$yhat.test.mean) <- rownames(x.test)
  names(res$prob.test.mean) <- rownames(x.test)

  res$nkeeptreedraws <- nkeeptreedraws
  res$ntree <- ntree
  res$call <- match.call()
  res$nobs <- n
  res$npred <- p
  class(res) <- "pmBART_fit"
  res
}

.validate_predictors <- function(x, name, allow_zero_rows = FALSE) {
  if (is.data.frame(x)) {
    if (!all(vapply(x, is.numeric, logical(1)))) {
      stop(sprintf("`%s` must contain only numeric columns.", name), call. = FALSE)
    }
    x <- as.matrix(x)
  }
  if (!is.matrix(x) || !is.numeric(x)) {
    stop(sprintf("`%s` must be a numeric matrix or data frame.", name), call. = FALSE)
  }
  if (ncol(x) < 1L) {
    stop(sprintf("`%s` must contain at least one predictor.", name), call. = FALSE)
  }
  if (!allow_zero_rows && nrow(x) < 2L) {
    stop(sprintf("`%s` must contain at least two observations.", name), call. = FALSE)
  }
  if (anyNA(x) || any(!is.finite(x))) {
    stop(sprintf("`%s` cannot contain missing or infinite values.", name), call. = FALSE)
  }
  storage.mode(x) <- "double"
  x
}

.check_scalar <- function(x, name, lower = -Inf, upper = Inf,
                          open_lower = FALSE) {
  valid_lower <- if (open_lower) x > lower else x >= lower
  if (length(x) != 1L || is.na(x) || !is.numeric(x) || !is.finite(x) ||
      !valid_lower || x > upper) {
    interval <- if (open_lower) "(" else "["
    stop(
      sprintf("`%s` must be a finite number in %s%s, %s].",
              name, interval, lower, upper),
      call. = FALSE
    )
  }
  invisible(x)
}

.check_positive_integer <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !is.numeric(x) ||
      !is.finite(x) || x < 1 || x != as.integer(x)) {
    stop(sprintf("`%s` must be a positive whole number.", name), call. = FALSE)
  }
  invisible(x)
}

.check_nonnegative_integer <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !is.numeric(x) ||
      !is.finite(x) || x < 0 || x != as.integer(x)) {
    stop(sprintf("`%s` must be a non-negative whole number.", name), call. = FALSE)
  }
  invisible(x)
}

.check_keep <- function(x, name, ndpost) {
  .check_nonnegative_integer(x, name)
  if (x > ndpost || (x > 0 && ndpost %% x != 0L)) {
    stop(
      sprintf("`%s` must be zero or a divisor of `ndpost` no larger than it.", name),
      call. = FALSE
    )
  }
  invisible(x)
}
