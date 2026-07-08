# Apply frozen quantile calibration to new samples

Maps test-sample feature values to the training distribution using the
frozen empirical-CDF grid from
[`fit_frozen_quantile`](https://kstawiski.github.io/OmicSelector/reference/fit_frozen_quantile.md).
Each test value is located in the training CDF, and the corresponding
training reference value is returned. This aligns the test distribution
to the training distribution in a monotone, feature-wise manner.

## Usage

``` r
apply_frozen_quantile(fit, x_test, missing_features = c("error", "skip"))
```

## Arguments

- fit:

  A `frozen_quantile_fit` object from
  [`fit_frozen_quantile`](https://kstawiski.github.io/OmicSelector/reference/fit_frozen_quantile.md).

- x_test:

  Numeric matrix (samples \\\times\\ features). All values must be
  non-negative.

- missing_features:

  `"error"` (default, fail-closed when training features are absent) or
  `"skip"` (process only shared features and return a reduced-column
  matrix).

## Value

Numeric matrix with columns for the shared features, values on the
training quantile scale.
