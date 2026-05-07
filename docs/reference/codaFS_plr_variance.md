# Pairwise Log-Ratio Variance Feature Filter

Scores pairwise log-ratios by combining their across-sample variance
with their association to the binary outcome, then aggregates the
highest-scoring ratios back to individual features by weighted
frequency.

## Usage

``` r
codaFS_plr_variance(
  X,
  y,
  n_features = NULL,
  positive = NULL,
  input_scale = c("log", "raw"),
  top_pairs = NULL,
  pair_multiplier = 5L,
  max_pair_features = 200L,
  prescreen_metric = c("clr_effect", "variance")
)
```

## Arguments

- X:

  Numeric matrix with samples in rows and features in columns.

- y:

  Binary outcome vector.

- n_features:

  Optional target panel size. Used only to derive the default number of
  pairwise log-ratios retained during aggregation.

- positive:

  Optional positive-class label for non-numeric outcomes.

- input_scale:

  Either \`"log"\` (default) or \`"raw"\`.

- top_pairs:

  Number of top pairwise log-ratios to aggregate. If \`NULL\`, this is
  set to \`pair_multiplier \* n_features\` or a conservative default.

- pair_multiplier:

  Multiplier used when \`top_pairs\` is \`NULL\`.

- max_pair_features:

  Maximum number of features retained before computing all pairwise
  log-ratios. Prescreening uses training-fold-only CLR effects.

- prescreen_metric:

  Either \`"clr_effect"\` or \`"variance"\`.

## Value

A list with per-feature \`scores\`, \`selected_features\`, and the
aggregated top pairwise log-ratio table.
