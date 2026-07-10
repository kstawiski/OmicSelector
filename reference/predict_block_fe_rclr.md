# Predict with a fitted block rCLR scorer

Predict with a fitted block rCLR scorer

## Usage

``` r
predict_block_fe_rclr(
  fit,
  expr,
  meta,
  cohort_col = "cohort",
  block_col = fit$block_col
)
```

## Arguments

- fit:

  Optional object from \`fit_block_fe_rclr()\`.

- expr:

  Numeric matrix, samples x features.

- meta:

  Data frame aligned to \`expr\`.

- cohort_col:

  Column in \`meta\` identifying cohorts; accepted for a common scorer
  signature but not used by this method.

- block_col:

  Column in \`meta\` identifying provenance blocks.

## Value

Numeric score vector.
