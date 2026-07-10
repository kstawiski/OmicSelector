# Fit the forced-single-sample entropic-OT (Sinkhorn) negative control

Fits the forced \\n = 1\\ entropic-OT negative control. The training
matrix is mapped to the per-sample robust CLR (over each row's own
positive support); the rows are split into a frozen CASE reference cloud
(`y_train == 1`) and a frozen CONTROL reference cloud (`y_train == 0`).
The frozen model stores only these two reference clouds (capped at
`hp$max_atoms` atoms each by a deterministic seeded subsample), the
frozen squared-cost scale, the entropic regularisation `hp$epsilon`, and
a single orientation scalar chosen so larger = more case-like. Scoring
is pure R and transports each specimen independently against these
frozen clouds.

## Usage

``` r
fit_sinkhorn_single(X_train, y_train, meta_train = NULL, hp = list())
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

  Optional list of hyperparameters. Allowed fields: `epsilon` (positive
  entropic-OT regularisation; default `0.1`), `min_features` (positive
  integer feature-overlap floor at scoring; default `3L`), `max_atoms`
  (positive integer cap on the atoms kept per reference cloud, seeded
  subsample; default `500L`), `eps` (positive finite numeric guard /
  scale floor; default `1e-8`), and `seed` (integer seed for the
  deterministic atom subsample; default `42L`).

## Value

Object of class `sinkhorn_single_model`: a list with `feature_universe`,
`case_atoms` / `control_atoms` (frozen rCLR reference clouds aligned to
`feature_universe`), `cost_scale`, `orientation`, `n_case` /
`n_control`, and `hp`.

## References

Cuturi M. (2013) Sinkhorn Distances. *NeurIPS* 26:2292-2300.

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
model <- fit_sinkhorn_single(X, y)
score_sinkhorn_single(model, X[1, , drop = FALSE])
} # }
```
