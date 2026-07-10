# Fit REO-singscore within-sample rank discriminator

Learns a frozen up/down feature signature from training data using
within-sample fractional ranks, then stores only that signature and
fixed hyperparameters for single-sample singscore-style scoring.
Training ranks are computed per sample with ascending ranks divided by
the training feature count; higher fractional rank means higher
abundance in that sample. A feature elevated in cases therefore has a
higher mean case fractional rank than control fractional rank, and is
selected by `mean_fracrank_cases - mean_fracrank_controls`. No test data
or test-time batch statistics are used.

## Usage

``` r
fit_reo_singscore(X_train, y_train, meta_train = NULL, hp = list())
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
  (default `min(25L, floor(ncol(X_train) / 4))`, at least 1) and
  `use_down` logical (default `TRUE`). These are resolved once during
  fitting and are not tuned inside `fit_reo_singscore()`.

## Value

A plain list of class `reo_singscore_model` containing `up_features`,
`down_features`, `feature_universe`, `use_down`, `k`, and the resolved
`hp`.

## References

Foroutan M, Bhuva DD, Lyu R, Horan K, Cursons J, Davis MJ. (2018) Single
sample scoring of molecular phenotypes. *BMC Bioinformatics* 19: 404.
PMID 30400809.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
X <- matrix(stats::rgamma(40 * 30, shape = 2), nrow = 40,
            dimnames = list(NULL, paste0("miR-", seq_len(30))))
y <- rep(c(0, 1), each = 20)
X[y == 1, 1:5] <- X[y == 1, 1:5] + 5
model <- fit_reo_singscore(X, y)
score <- score_reo_singscore(model, X)
} # }
```
