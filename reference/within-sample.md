# Within-Sample Normalization for Biomarker Panels

Functions for within-sample normalization that operate ONLY on a single
sample's feature values. Because these methods use no population-level
statistics — each sample is normalized independently using only its own
feature values — they cancel specific classes of per-sample nuisance
variation WITHOUT a reference cohort. Invariance scope is
method-specific (see Details). These methods are NOT a general
"batch-effect remover": gene-specific batch effects, non-linear platform
response, and measurement-specific bias are NOT cancelled.

This is the key insight for clinical deployment: a diagnostic test based
on within-sample normalization requires NO reference cohort. Draw blood,
measure k biomarkers, compute the relative pattern, classify.

## Details

Within-sample normalization addresses one component of the batch-effect
problem in circulating biomarker diagnostics: per-sample nuisance
variation (collection-site offset, loading amount, total-library
scaling). Traditional normalization methods (quantile, ComBat, z-score)
require a reference population, which is impractical for point-of-care
diagnostics. Within-sample methods sidestep the reference-cohort
requirement but each has a specific, narrower invariance scope:

- `ws_logratio` / `ws_ratio_image`: additive-shift-invariant in
  log-space (equivalently, invariant to a multiplicative per-sample
  factor applied UNIFORMLY across all features in raw space). NOT
  invariant to gene-specific shifts or non-linear platform response.

- `ws_zscore`: invariant to affine per-sample transforms (scale + shift)
  applied uniformly across features.

- `ws_rank`: invariant to any strictly monotonic per-sample transform of
  all features.

- `ws_minmax`: invariant to uniform per-sample affine transforms but
  sensitive to outlier features.

Available methods:

- `ws_minmax`: Scale each sample's features to \[0,1\]

- `ws_rank`: Replace values with within-sample ranks

- `ws_zscore`: Standardize across features within each sample

- `ws_logratio`: Compute all pairwise log-ratios (self-normalizing)

- `ws_ratio_image`: Create pairwise ratio matrix (for CNN input)
