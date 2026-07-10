# Predict biofluid-stable anchor rCLR scores

Predict biofluid-stable anchor rCLR scores

## Usage

``` r
predict_biofluid_anchor_rclr(
  fit,
  expr,
  meta,
  biofluid_col = fit$biofluid_col,
  feature_weights = NULL
)
```

## Arguments

- fit:

  Object from \`fit_biofluid_anchor_rclr()\`.

- expr:

  Numeric matrix, samples x features.

- meta:

  Data frame aligned to \`expr\`.

- biofluid_col:

  Column in \`meta\` with canonical or raw biofluid.

- feature_weights:

  Optional named numeric panel weights.
