# Fit kit-stratified rCLR centering sets.

Fit kit-stratified rCLR centering sets.

## Usage

``` r
fit_kit_stratified_rclr(
  expr_matrix,
  sample_meta,
  kit_label_col,
  biofluid_col = "biofluid",
  trim_upper = 0.1,
  trim_lower = 0.05,
  min_centering_size = 8L,
  min_samples_per_stratum = 5L,
  exclude_features = c("hsa-miR-451a", "hsa-miR-16-5p", "hsa-miR-486-5p",
    "hsa-miR-144-3p", "hsa-miR-223-3p")
)
```

## Arguments

- expr_matrix:

  Numeric matrix, samples x features, counts or non-negative abundance
  values.

- sample_meta:

  data.frame with one row per sample.

- kit_label_col:

  Column in \`sample_meta\` containing library-kit family.
  Missing/unknown kit values fall back to \`biofluid_col\`, then global.

- biofluid_col:

  Optional fallback column. Default "biofluid".

- trim_upper, trim_lower:

  Feature trimming fractions used to choose centering features from the
  training stratum.

- min_centering_size:

  Minimum centering set size before global fallback.

- min_samples_per_stratum:

  Minimum training samples needed for a stratum-specific centering set.

- exclude_features:

  Features excluded from denominator selection.

## Value

A fit object consumed by \`score_kit_stratified_rclr()\`.
