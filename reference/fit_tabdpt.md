# Fit the TabDPT in-context discriminator (freeze the training context)

Freezes the training context for TabDPT in-context classification. The
training matrix is mapped to the per-sample robust CLR over the frozen
feature universe (`colnames(X_train)`); if the context exceeds
`hp$max_context` rows it is reduced to a deterministic, seeded,
class-stratified subsample. A TabDPT classifier (pinned token-free
checkpoint) is constructed and `fit()` on that frozen rCLR context
(which, for TabDPT, simply loads the context into the model and fits the
preprocessing scaler / optional PCA – no gradient training). The frozen
R-side state needed to reproduce scores (the rCLR context, labels,
feature universe, class order, device/config, hp) is stored on the
model; the live python classifier is also kept for scoring but is NEVER
part of the model digest.

## Usage

``` r
fit_tabdpt(X_train, y_train, meta_train = NULL, hp = list())
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

  Optional list of hyperparameters. Allowed fields: `device` (`"cpu"`
  (default), `"cuda"`, or `"auto"`; `"cuda"` / `"auto"` fall back to CPU
  when no GPU is available), `temperature` (TabDPT softmax temperature
  at predict, positive finite, default `1.0`), `normalizer` (TabDPT
  input scaler, one of `"standard"` (default), `"minmax"`, `"robust"`,
  `"power"`, `"quantile-uniform"`, `"quantile-normal"`, `"log1p"`),
  `max_context` (row cap on the frozen context, integer \\\ge 2\\,
  default `4096L`; larger contexts are seeded class-stratified
  subsampled), `min_features` (feature-overlap floor at scoring,
  positive integer, default `3L`), `seed` (integer; default `42L`), and
  `score_batch` (benchmark-only batched scoring flag, logical, default
  `FALSE`).

## Value

Object of class `tabdpt_model`: a list with `context_rclr` (frozen rCLR
context matrix), `context_y`, `feature_universe`, `classes`, `case_col`,
`n_context`, `device` (resolved), `checkpoint` (pinned checkpoint
version), `hp`, and the live python classifier `clf` (excluded from the
digest).

## References

Ma J, Thomas V, Hosseinzadeh R, et al. (2024) TabDPT: Scaling Tabular
Foundation Models. arXiv:2410.18164.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
n <- 120; p <- 20; k <- 6
L <- matrix(stats::rnorm(n * p, 4, 0.5), nrow = n,
            dimnames = list(NULL, paste0("miR-", seq_len(p))))
y <- rep(c(0, 1), each = n / 2)
L[y == 1, seq_len(k)] <- L[y == 1, seq_len(k)] + 1.2
X <- exp(L)
model <- fit_tabdpt(X, y)
score_tabdpt(model, X[1, , drop = FALSE])
} # }
```
