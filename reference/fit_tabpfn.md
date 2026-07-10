# Fit the TabPFN-v2 in-context discriminator (freeze the training context)

Freezes the training context for TabPFN-v2 in-context classification.
The training matrix is mapped to the per-sample robust CLR over the
frozen feature universe (`colnames(X_train)`); if the context exceeds
`hp$max_context` rows it is reduced to a deterministic, seeded,
class-stratified subsample. A TabPFN-v2 (ungated `ModelVersion.V2`)
classifier is constructed and `fit()` on that frozen rCLR context
(which, for TabPFN, simply loads the context into the model – no
gradient training). The frozen R-side state needed to reproduce scores
(the rCLR context, labels, feature universe, class order, device/config,
hp) is stored on the model; the live python classifier is also kept for
scoring but is NEVER part of the model digest. TabPFN-v2 is applied to
the full rCLR feature universe; cohorts with \\\>500\\ features,
including approximately 2540-2565-feature high-plex miRNA panels, exceed
TabPFN's nominal 500-feature ceiling, so its documented
`ignore_pretraining_limits` override is enabled. Under the pinned
`ModelVersion.V2` preprocessor path (tabpfn 8.0.7) the per-estimator
feature-subsampling threshold is 1,000,000 (the 500-feature
random-subspace reduction is the separate V2.5 preset, which is not used
here), so TabPFN-v2 ingests the full feature set directly and is run
outside its nominal \\\le 500\\-feature pretraining regime. No
per-estimator feature subsampling is introduced (the ensemble subsample
indices are all empty), so scoring stays deterministic and exactly
row-equivariant. The override also lifts the CPU \\\>1000\\-sample guard
(a compute guard, not a numeric change). This matches how the sibling
in-context foundation models TabICL and TabDPT are run on the full
representation and keeps the comparison apples-to-apples with the rCLR
baseline; a leakage-safe in-regime top-500 feature reduction is
evaluated as a companion sensitivity analysis (the prespecified
out-of-regime control).

## Usage

``` r
fit_tabpfn(X_train, y_train, meta_train = NULL, hp = list())
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
  when no GPU is available), `n_estimators` (TabPFN ensemble size,
  positive integer, default `8L`), `max_context` (row cap on the frozen
  context, integer \\\ge 2\\, default `4096L`; larger contexts are
  seeded class-stratified subsampled), `min_features` (feature-overlap
  floor at scoring, positive integer, default `3L`), `seed` (integer;
  default `42L`), and `score_batch` (benchmark-only batched scoring
  flag, logical, default `FALSE`).

## Value

Object of class `tabpfn_model`: a list with `context_rclr` (frozen rCLR
context matrix), `context_y`, `feature_universe`, `classes`, `case_col`,
`n_context`, `device` (resolved), `model_version` (`"v2"`),
`ignore_pretraining_limits` (`TRUE`), `hp`, and the live python
classifier `clf` (excluded from the digest).

## References

Hollmann N, et al. (2025) Accurate predictions on small data with a
tabular foundation model. *Nature* 637:319-326.

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
model <- fit_tabpfn(X, y)
score_tabpfn(model, X[1, , drop = FALSE])
} # }
```
