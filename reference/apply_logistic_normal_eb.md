# Apply frozen logistic-normal EB shrinkage to new samples

Shrinks CLR-transformed test samples toward the per-feature training
prior using the model fitted by
[`fit_logistic_normal_eb`](https://kstawiski.github.io/OmicSelector/reference/fit_logistic_normal_eb.md).

Missing test features are handled by `missing_features`: `"error"`
(default, fail-closed) or `"fill_with_prior"` (pads missing features
with `NA` in the test matrix, then the posterior collapses fully to the
prior mean for those features).

## Usage

``` r
apply_logistic_normal_eb(
  fit,
  x_test,
  missing_features = c("error", "fill_with_prior")
)
```

## Arguments

- fit:

  A `logistic_normal_eb_fit` object from
  [`fit_logistic_normal_eb`](https://kstawiski.github.io/OmicSelector/reference/fit_logistic_normal_eb.md).

- x_test:

  Numeric matrix (samples \\\times\\ features). All values must be
  non-negative.

- missing_features:

  `"error"` (default) or `"fill_with_prior"`.

## Value

Numeric matrix (same row/col structure as ordered training features) of
posterior-shrunk CLR coordinates.
