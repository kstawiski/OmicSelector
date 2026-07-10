# Fit the CoDaCoRe stagewise log-ratio-balance discriminator

Trains the CoDaCoRe continuous relaxation stagewise (boosting) via
reticulate-python torch (CPU by default for exact reproducibility),
DISCRETIZES each balance to two disjoint feature sets, FREEZES them,
fits a base-R IRLS ridge logistic head over the discrete balances,
EXPORTS the discrete sets + head to R, and DISCARDS the python module.
The fitted model holds no external pointer.

## Usage

``` r
fit_coda_codacore(X_train, y_train, meta_train = NULL, hp = list())
```

## Arguments

- X_train:

  Numeric matrix (samples \\\times\\ features) of non-negative
  abundances with unique, non-empty feature names (colnames).

- y_train:

  Numeric / integer 0/1 labels (1 = case), length `nrow(X_train)`, with
  at least one case and one control.

- meta_train:

  Optional per-sample metadata. Accepted for interface uniformity and
  otherwise ignored.

- hp:

  Optional list of hyperparameters. Allowed fields: `max_balances` (max
  stagewise balances, positive integer; default `3L`), `epochs` (Adam
  epochs per stage, positive integer; default `200L`), `lr` (positive
  Adam learning rate; default `1e-2`), `tau` (positive softmax
  relaxation temperature; default `1.0`), `ridge` (non-negative L2 ridge
  on the logistic head; default `1e-4`), `min_features` (feature-overlap
  floor at scoring, positive integer; default `3L`), `device` (`"cpu"`
  (default), `"cuda"`, or `"auto"`; `"cuda"`/`"auto"` fall back to CPU
  with no GPU – but note the default `"cpu"` is what guarantees an
  EXACTLY reproducible discretization), and `seed` (integer; default
  `42L`).

## Value

Object of class `coda_codacore_model`: a list with `feature_universe`,
`balances` (a list of `list(A, B)` 1-based universe-position index sets
per frozen discrete balance), `intercept`, `weights` (the frozen
logistic head over the balances), `n_balances`, `device` (resolved),
`seed`, and `hp`. No python pointer.

## References

Gordon-Rodriguez E, Quinn TP, Cunningham JP. (2021) Learning sparse
log-ratios for high-throughput sequencing data. *Bioinformatics*
38(1):157-163.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
n <- 120; p <- 16; k <- 3
X <- matrix(stats::rgamma(n * p, shape = 30, rate = 2), nrow = n,
            dimnames = list(NULL, paste0("f", sprintf("%02d", seq_len(p)))))
y <- rep(c(0, 1), each = n / 2)
X[y == 1, 1:k] <- X[y == 1, 1:k] * 1.8
X[y == 1, (p - k + 1):p] <- X[y == 1, (p - k + 1):p] / 1.8
model <- fit_coda_codacore(X, y)
score_coda_codacore(model, X[1, , drop = FALSE])
} # }
```
