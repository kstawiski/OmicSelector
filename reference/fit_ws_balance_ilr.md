# Fit within-cohort ILR balance discriminator

Wraps
[`ws_balance_ilr`](https://kstawiski.github.io/OmicSelector/reference/ws_balance_ilr.md)
as a frozen single-sample discriminator. The raw score is the row sum of
valid sequential-binary-partition ILR balances computed from each
specimen's own feature values. Fitting learns only a one-bit orientation
from the training labels: if the mean raw score is at least as high in
cases as controls, the stored sign is `+1`; otherwise it is `-1`. No
test data or test-time batch statistic is used.

## Usage

``` r
fit_ws_balance_ilr(X_train, y_train, meta_train = NULL, hp = list())
```

## Arguments

- X_train:

  Numeric matrix (samples \\\times\\ features) of non-negative abundance
  values. Columns must be uniquely named feature ids.

- y_train:

  Integer/numeric 0/1 labels aligned to `X_train`; 1 is case/disease and
  0 is control.

- meta_train:

  Optional per-sample metadata. Accepted for interface uniformity and
  ignored by this method.

- hp:

  List of frozen hyperparameters: `sbp` sequential binary partition
  (default
  [`ws_default_sbp()`](https://kstawiski.github.io/OmicSelector/reference/ws_default_sbp.md)),
  `pseudocount` (default `NULL`, using the primitive's per-sample
  default), `aggregate` (`"gmean"` or `"trimmed_gmean"`, default
  `"gmean"`), `min_balance_coverage` (default `0.8`), and
  `min_sbp_features` canonical SBP feature floor (default `3L`).

## Value

A plain list of class `ws_balance_ilr_model` containing the frozen
`sbp`, `pseudocount`, `aggregate`, `min_balance_coverage`,
`min_sbp_features`, orientation `sign`, `sbp_canonical_features`,
`feature_universe`, and resolved `hp`.

## References

Egozcue J. J., Pawlowsky-Glahn V. (2005) Groups of parts and their
balances in compositional data analysis. *Mathematical Geology* 37(7):
795-828.

## Examples

``` r
if (FALSE) { # \dontrun{
sbp <- ws_default_sbp()
features <- unique(unlist(lapply(sbp, function(b) {
  c(b$numerator, b$denominator)
}), use.names = FALSE))
features <- setdiff(features, c("ALL_NON_RBC", "ALL_NON_PLT",
                                "TOP_K_BY_ABUNDANCE", "TAIL_BY_ABUNDANCE"))
X <- matrix(stats::rgamma(40 * length(features), shape = 2), nrow = 40,
            dimnames = list(NULL, features))
y <- rep(c(0, 1), each = 20)
model <- fit_ws_balance_ilr(X, y)
score <- score_ws_balance_ilr(model, X)
} # }
```
