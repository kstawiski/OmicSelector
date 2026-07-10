# Fit the truncated path-signature single-sample discriminator

Maps the TRAINING matrix to the per-sample robust CLR over the frozen
feature universe (`colnames(X_train)`), computes each training
specimen's closed-form truncated path signature (time-augmented 2D path,
level `L`, pure base R, ROW-BY-ROW), freezes a per-coordinate
signature-coefficient standardization, and fits a FROZEN ridge-logistic
head over the standardized coefficients. If `hp$lambda` is `NULL` the
ridge penalty is selected by TRAINING-ONLY deterministic
class-stratified k-fold CV (no RNG). The fitted model holds NO python
dependency and NO external pointer.

## Usage

``` r
fit_sig_path(X_train, y_train, meta_train = NULL, hp = list())
```

## Arguments

- X_train:

  Numeric matrix (samples \\\times\\ features) of non-negative
  abundances with unique, non-empty feature names (colnames). At least 2
  features are required (a path needs at least one segment).

- y_train:

  Numeric / integer 0/1 labels (1 = case), length `nrow(X_train)`, with
  at least one case and one control.

- meta_train:

  Optional per-sample metadata. Accepted for interface uniformity and
  otherwise ignored.

- hp:

  Optional list of hyperparameters. Allowed fields: `L` (signature
  truncation level, one of `1L`, `2L`, `3L`, `4L`; default `3L`;
  signature dimension \\D = \sum\_{k=1}^{L} 2^k\\, so `14` at `L = 3`),
  `lambda` (ridge penalty for the logistic head; `NULL` (default)
  selects it by training-only deterministic k-fold CV, otherwise a
  non-negative finite number), `min_features` (feature-overlap floor at
  scoring, positive integer; default `3L`), and `seed` (integer; default
  `42L`). NO `device` hp — the method is pure base R.

## Value

Object of class `sig_path_model`: a list with `feature_universe`, `L`,
`D` (signature dimension), `coef_mean`, `coef_sd` (frozen
standardization), `intercept`, `weights` (frozen linear head), `lambda`
(resolved ridge penalty), `seed`, and `hp`. No python dependency, no
external pointer.

## References

Chevyrev I, Kormilitzin A. (2016) A primer on the signature method in
machine learning. *arXiv:1603.03788*.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
n <- 80; p <- 16; k <- 6
L <- matrix(stats::rnorm(n * p, 4, 0.5), nrow = n,
            dimnames = list(NULL, paste0("miR-", seq_len(p))))
y <- rep(c(0, 1), each = n / 2)
L[y == 1, seq_len(k)] <- L[y == 1, seq_len(k)] + 1.2
X <- exp(L)
model <- fit_sig_path(X, y)
score_sig_path(model, X[1, , drop = FALSE])
} # }
```
