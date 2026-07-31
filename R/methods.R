#' @export
print.pmBART_fit <- function(x, ...) {
  cat("Probit monotone BART fit\n")
  cat("  Observations:", x$nobs, "\n")
  cat("  Predictors:  ", x$npred, "\n")
  cat("  Trees:       ", x$ntree, "\n")
  cat("  Tree draws:  ", x$nkeeptreedraws, "\n")
  if (ncol(x$prob.test) > 0L) {
    cat("  Test rows:   ", ncol(x$prob.test), "\n")
  }
  invisible(x)
}
