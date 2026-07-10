# Predict with a fitted cohort-z scorer

Predict with a fitted cohort-z scorer

## Usage

``` r
predict_cohort_z_then_pool_score(
  fit,
  expr,
  meta,
  cohort_col = fit$cohort_col,
  block_col = "provenance_block",
  allow_new_group_moments = TRUE
)
```

## Arguments

- fit:

  Optional object from \`fit_cohort_z_then_pool_score()\`.

- expr:

  Numeric matrix, samples x features.

- meta:

  Data frame aligned to \`expr\`.

- cohort_col:

  Column in \`meta\` identifying cohorts.

- block_col:

  Column in \`meta\` identifying provenance blocks; accepted for a
  common scorer signature but not used by this method.

- allow_new_group_moments:

  Logical; if TRUE, new cohorts are centered on their own unlabeled
  samples.

## Value

Numeric score vector.
