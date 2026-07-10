# Score kit-orthogonal ILR-like residualized panel.

Score kit-orthogonal ILR-like residualized panel.

## Usage

``` r
score_kit_orthogonal_ilr(
  expr_matrix,
  sample_meta,
  panel_features,
  kit_label_col,
  biofluid_col = "biofluid",
  fit = NULL,
  feature_weights = NULL,
  pseudocount = NULL,
  ...
)
```

## Arguments

- expr_matrix:

  Numeric matrix, samples x features.

- sample_meta:

  data.frame with one row per sample and \`kit_label_col\`.

- panel_features:

  Fold-specific panel features, typically top 20.

- kit_label_col:

  Column containing kit family. NA, "unknown", "NR", and non-NGS
  sentinels fall back to biofluid, then global centering.

- biofluid_col:

  Optional fallback stratum column.

- fit:

  Optional object from \`fit_kit_orthogonal_ilr()\`.

- feature_weights:

  Optional named numeric signs/weights for panel features; default all
  +1. This lets training folds orient features so higher scores are more
  case-like.

- pseudocount:

  Optional additive pseudocount before log transform.

- ...:

  Additional arguments passed to the fitting helper when \`fit\` is not
  supplied.

## Value

Numeric score vector of length nrow(expr_matrix).
