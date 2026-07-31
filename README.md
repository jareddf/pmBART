# pmBART

`pmBART` fits probit monotone BART models for binary outcomes. It constrains
every predictor to have a non-decreasing effect on the latent response.

This package is intended for installation from source rather than submission to
CRAN.

## Installation

Install the package from its local directory:

```r
install.packages(
  "/path/to/pmBART",
  repos = NULL,
  type = "source"
)
```

For development, `devtools::load_all()` is convenient.

## Quick start

```r
library(pmBART)

set.seed(42)
x <- cbind(
  exposure = runif(100),
  score = runif(100)
)
probability <- pnorm(-1 + 1.5 * x[, "exposure"] + x[, "score"])
y <- rbinom(100, size = 1, prob = probability)

fit <- probit_monbart(
  x.train = x[1:80, ],
  y.train = y[1:80],
  x.test = x[81:100, ]
)

fit
fit$prob.train.mean
fit$prob.test.mean
```

The two `prob.*.mean` vectors are the simplest outputs for fitted values and
test-set predictions. The corresponding `prob.train` and `prob.test` matrices
contain retained posterior draws and can be used to calculate uncertainty
intervals:

```r
test_intervals <- apply(
  fit$prob.test,
  2,
  quantile,
  probs = c(0.025, 0.5, 0.975)
)
```

## Important modeling assumptions

- All predictors must be numeric and cannot contain missing values.
- The response must contain both `0` and `1`.
- Every predictor is constrained to have a non-decreasing effect. To encode a
  decreasing effect, reverse that predictor before fitting (for example, use
  `-x`).
- New prediction points are supplied with `x.test` when the model is fitted.
- Use `?probit_monbart` for the complete argument and return-value reference.
