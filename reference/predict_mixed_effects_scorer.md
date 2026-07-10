# Predict a mixed-effects panel score without reading test labels

Predict a mixed-effects panel score without reading test labels

## Usage

``` r
predict_mixed_effects_scorer(
  fit,
  expr,
  meta,
  cohort_col = fit$cohort_col,
  block_col = fit$block_col,
  min_weighted_features = 3L
)
```

## Arguments

- fit:

  Object from \`fit_mixed_effects_scorer()\`.

- expr:

  Numeric matrix, test samples x features.

- meta:

  Data frame aligned to \`expr\`; labels are not used.

- cohort_col:

  Test cohort column.

- block_col:

  Test provenance-block column.

- min_weighted_features:

  Minimum finite weights required for scoring.

## Value

Numeric score vector, one value per sample.
