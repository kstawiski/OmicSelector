# Predict biofluid-residualized ALR scores

Predict biofluid-residualized ALR scores

## Usage

``` r
predict_biofluid_residualized_alr(
  fit,
  expr,
  meta,
  biofluid_col = fit$biofluid_col,
  feature_weights = NULL
)
```

## Arguments

- fit:

  Object from \`fit_biofluid_residualized_alr()\`.

- expr:

  Numeric matrix, samples x features.

- meta:

  Data frame aligned to \`expr\`.

- biofluid_col:

  Column identifying biofluid class.

- feature_weights:

  Optional named numeric panel weights.
