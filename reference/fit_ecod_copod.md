# Fit the ECOD + COPOD novelty discriminator

Fits the bundled ECOD + COPOD novelty discriminator. The training matrix
is mapped to the per-sample robust CLR (over each row's own positive
support), and the in-distribution reference is taken to be the training
CONTROL rows (`y_train == 0`). pyod's `ECOD` and `COPOD` are fitted on
that control rCLR (leakage-safe: training only) and used to (a) source /
verify the per-feature skewness signs and (b) provide the
`decision_scores_` anchor that the frozen pure-R reimplementation must
reproduce. The frozen model stores only R-side state – the per-feature
control reference columns, the skew signs, the per-detector
control-aggregate mean / sd, and the orientation scalar – so scoring is
pure R and never calls Python.

## Usage

``` r
fit_ecod_copod(X_train, y_train, meta_train = NULL, hp = list())
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

  Optional list of hyperparameters. Allowed fields: `min_features`
  (positive integer feature-overlap floor at scoring; default `3L`),
  `eps` (positive finite floor for the standardization sd and a numeric
  guard; default `1e-8`), and `verify_pyod` (single logical; when `TRUE`
  (default) the frozen pure-R reimplementation is checked against pyod's
  `decision_scores_` on the control rows at fit and the fit fails if the
  maximum difference exceeds `1e-6`).

## Value

Object of class `ecod_copod_model`: a list with `feature_universe`,
`ref_cols` (frozen per-feature control reference columns), `skew_sign`,
`score_mean` / `score_sd` (per-detector control-aggregate
standardization constants), `orientation`, `n_control`, and `hp`.

## References

Li Z, Zhao Y, Hu X, Botta N, Ionescu C, Chen G. (2022) ECOD. *IEEE
TKDE*.

Li Z, Zhao Y, Botta N, Ionescu C, Hu X. (2020) COPOD. *ICDM*.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
n <- 120; p <- 30; k <- 8
L <- matrix(stats::rnorm(n * p, 4, 0.5), nrow = n,
            dimnames = list(NULL, paste0("miR-", seq_len(p))))
y <- rep(c(0, 1), each = n / 2)
L[y == 1, seq_len(k)] <- L[y == 1, seq_len(k)] + 1.2
X <- exp(L)
model <- fit_ecod_copod(X, y)
score_ecod_copod(model, X[1, , drop = FALSE])
} # }
```
