# Fit REO-UCell within-sample rank discriminator

Learns a frozen up/down feature signature from training data using
within-sample ranks, then stores only that signature and fixed
hyperparameters for single-sample UCell-style scoring. Training ranks
are computed per sample with rank 1 assigned to the most abundant
feature. A feature elevated in cases therefore has a lower mean case
rank than control rank, and is selected by
`mean_rank_controls - mean_rank_cases`. No test data or test-time batch
statistics are used.

## Usage

``` r
fit_reo_ucell(X_train, y_train, meta_train = NULL, hp = list())
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

  List of frozen hyperparameters: `k` signature size per direction
  (default `min(25L, floor(ncol(X_train) / 4))`, at least 1), `maxRank`
  UCell rank cap (default `min(1500L, ncol(X_train))`), and `use_down`
  logical (default `TRUE`). These are resolved once during fitting and
  are not tuned inside `fit_reo_ucell()`.

## Value

A plain list of class `reo_ucell_model` containing `up_features`,
`down_features`, `maxRank`, `feature_universe`, `use_down`, `k`, and the
resolved `hp`.

## References

Andreatta M, Carmona SJ. (2021) UCell: Robust and scalable single-cell
gene signature scoring. *Computational and Structural Biotechnology
Journal* 19: 3796-3798.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
X <- matrix(stats::rgamma(40 * 30, shape = 2), nrow = 40,
            dimnames = list(NULL, paste0("miR-", seq_len(30))))
y <- rep(c(0, 1), each = 20)
X[y == 1, 1:5] <- X[y == 1, 1:5] + 5
model <- fit_reo_ucell(X, y)
score <- score_reo_ucell(model, X)
} # }
```
