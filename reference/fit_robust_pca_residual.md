# Fit a MAD-scaled SVD residual filter for batch correction

Fits a frozen MAD-scaled PCA residual filter. This is NOT the Candès
RPCA (ALM sparse-low-rank decomposition) — it is classical SVD on a
robustified pre-processing step (per-feature median centering + MAD
scaling), which is more resistant to outlier samples than vanilla
CLR-then-PCA but is not the L+S sparse-low-rank framework.

Training pipeline: per-feature median + MAD scaling, then classical SVD;
retain top-`n_components` right singular vectors as a frozen PC basis.

## Usage

``` r
fit_robust_pca_residual(x_train, n_components = 5L, pseudocount = 0.5)
```

## Arguments

- x_train:

  Numeric matrix (samples \\\times\\ features). Must have column names.
  All values must be non-negative.

- n_components:

  Integer — number of principal components to retain (and remove).
  Default 5.

- pseudocount:

  Numeric \\\> 0\\. Added before log transform. Default 0.5.

## Value

Object of class `robust_pca_residual_fit` with components:
`feature_names`, `feat_median`, `feat_mad`, `v` (right singular vectors,
features \\\times\\ `n_components`), `n_components`, `pseudocount`,
`method`.

## References

Hubert M, Rousseeuw PJ, Vanden Branden K. (2005) ROBPCA. *Technometrics*
47(1): 64–79.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
x_train <- matrix(abs(rnorm(100 * 30, mean = 100, sd = 30)), nrow = 100,
                  dimnames = list(NULL, paste0("miR-", seq_len(30))))
fit <- fit_robust_pca_residual(x_train, n_components = 3L)
} # }
```
