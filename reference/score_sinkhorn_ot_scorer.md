# Score projected samples by summing panel features.

Score projected samples by summing panel features.

## Usage

``` r
score_sinkhorn_ot_scorer(
  X_test,
  fit,
  panel_features = fit$panel_features,
  feature_weights = fit$feature_weights,
  meta_test = NULL
)
```

## Arguments

- X_test:

  Samples x features matrix in the fitted feature space.

- fit:

  Object from \`fit_sinkhorn_ot_scorer()\`.

- panel_features:

  Panel features to sum after projection.

- feature_weights:

  Optional named numeric panel weights.

- meta_test:

  Optional metadata retained for a common scorer signature.
