# Apply compositional Mahalanobis distance to new samples

Projects new samples into the same log-ratio space and computes their
squared Mahalanobis distance using the frozen center and
inverse-covariance from
[`fit_compositional_mahalanobis`](https://kstawiski.github.io/OmicSelector/reference/fit_compositional_mahalanobis.md).

## Usage

``` r
apply_compositional_mahalanobis(fit, x_test, return_pvalue = TRUE)
```

## Arguments

- fit:

  A `compositional_mahalanobis_fit` object.

- x_test:

  Numeric matrix (samples \\\times\\ features). All values must be
  non-negative. All training features must be present.

- return_pvalue:

  Logical. If `TRUE` (default), also returns chi-square reference
  p-values.

## Value

List with `distance_sq`, `distance`, and (if `return_pvalue`)
`pvalue_chisq`.
