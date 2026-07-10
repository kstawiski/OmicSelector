# Fit a Group-DRO kit x biofluid logistic scorer.

Fit a Group-DRO kit x biofluid logistic scorer.

## Usage

``` r
fit_group_dro_scorer(
  expr_matrix,
  y,
  sample_meta,
  kit_col = "kit_family",
  biofluid_col = "biofluid",
  panel_features = NULL,
  panel_size = 20L,
  eta_grid = c(0.01, 0.05, 0.1),
  l2_grid = c(0.001, 0.01, 0.1),
  inner_folds = 3L,
  learning_rate = 0.05,
  epochs = 200L,
  tune = TRUE,
  seed = 42L,
  ...
)
```

## Arguments

- expr_matrix:

  Numeric training matrix, samples x miRNA features.

- y:

  Binary 0/1 training labels.

- sample_meta:

  Metadata aligned to rows of \`expr_matrix\`.

- kit_col, biofluid_col:

  Metadata columns used for train-only group construction. Groups are
  \`kit x biofluid\`.

- panel_features:

  Optional training-selected feature panel. If omitted, a train-only CLR
  panel is selected inside this function.

- panel_size:

  Number of features selected when \`panel_features\` is NULL.

- eta_grid, l2_grid:

  Hyperparameter grids for train-only nested CV.

- inner_folds:

  Number of inner folds for hyperparameter tuning.

- learning_rate:

  Gradient-descent learning rate.

- epochs:

  Number of optimization epochs.

- tune:

  Logical. If TRUE, tune \`eta\` and \`l2\` on training data.

- seed:

  Integer random seed.

- ...:

  Reserved for scorer-interface compatibility.

## Value

A frozen Group-DRO fit object consumed by \`score_group_dro_scorer()\`.
