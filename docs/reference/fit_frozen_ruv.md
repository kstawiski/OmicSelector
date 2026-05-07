# Fit a frozen RUV factor model from training data

Estimates latent factors of unwanted variation from empirical
negative-control features on training data, then freezes the factor
loadings for deployment. This is the frozen-mode extension of RUVSeq
(Risso et al., Nat Biotechnol 2014): for a new sample \\x\\, the
unwanted variation component is estimated via the frozen pseudo-inverse
of training latent factors, then subtracted.

If fewer than 5 named control features are found in `colnames(X_train)`,
falls back to the top-30% lowest-variance features as proxy controls,
with a warning.

## Usage

``` r
fit_frozen_ruv(
  X_train,
  control_features = NULL,
  k = 2:6,
  select = c("BIC", "fixed", "explicit")
)
```

## Arguments

- X_train:

  Numeric matrix (samples \\\times\\ features). Training cohort.

- control_features:

  Character vector of feature names expected to be non-cancer-relevant
  (empirical negative controls). Default = a curated set of
  stably-expressed serum/plasma miRNAs.

- k:

  Integer or integer vector (range). Number of unwanted factors to
  estimate. Default `2:6` — BIC selects within this range.

- select:

  Method for selecting `k`. `"BIC"` (default), `"fixed"` (uses the
  smallest value in `k`), or `"explicit"` (if `select` is a single
  numeric).

## Value

Object of class `frozen_ruv_fit` with components: `feat_mean`,
`full_loadings`, `latent_train_inv`, `k_selected`, `bics`, `n_train`,
`p_full`, `n_controls`.

## References

Risso D, Ngai J, Speed TP, Dudoit S. (2014) Normalization of RNA-seq
data using factor analysis of control genes or samples. *Nature
Biotechnology* 32(9): 896–902.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
X <- matrix(rnorm(200 * 50), nrow = 200,
            dimnames = list(NULL, paste0("mir-", seq_len(50))))
fit <- fit_frozen_ruv(X, k = 2:4)
fit$k_selected
} # }
```
