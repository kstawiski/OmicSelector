# Apply MAD-scaled SVD residual filter to new samples

Projects test samples into the frozen PC space estimated by
[`fit_robust_pca_residual`](https://kstawiski.github.io/OmicSelector/reference/fit_robust_pca_residual.md),
then returns the residual (test sample minus reconstruction from
top-`n_components` PCs) as a low-rank-removed representation.

All training features must be present in `x_test`.

## Usage

``` r
apply_robust_pca_residual(fit, x_test, return_reconstructed = FALSE)
```

## Arguments

- fit:

  A `robust_pca_residual_fit` object from
  [`fit_robust_pca_residual`](https://kstawiski.github.io/OmicSelector/reference/fit_robust_pca_residual.md).

- x_test:

  Numeric matrix (samples \\\times\\ features). All values must be
  non-negative.

- return_reconstructed:

  Logical. If `TRUE`, returns a list with components `scores`,
  `reconstructed`, and `residual`. If `FALSE` (default), returns only
  the residual matrix.

## Value

The residual matrix (samples \\\times\\ features) on the MAD-scaled log
scale, or a list when `return_reconstructed = TRUE`.
