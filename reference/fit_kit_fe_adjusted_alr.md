# Fit kit fixed-effect adjusted ALR.

Per feature, log-CPM values are adjusted by subtracting the training
stratum mean and adding the global training mean. ALR pivots are
selected from features with the lowest between-kit variance. When kit is
unknown, biofluid and then global strata are used.

## Usage

``` r
fit_kit_fe_adjusted_alr(
  expr_matrix,
  sample_meta,
  kit_label_col,
  biofluid_col = "biofluid",
  n_pivots = 6L,
  min_samples_per_stratum = 5L,
  pseudocount = 0.5,
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

- n_pivots:

  Number of low-variance pivot features to select.

- min_samples_per_stratum:

  Minimum training samples needed for a stratum-specific centering set.

- pseudocount:

  Additive pseudocount before log-CPM transformation.

- exclude_features:

  Features excluded from denominator selection.

## Value

Fit object consumed by \`score_kit_fe_adjusted_alr()\`.
