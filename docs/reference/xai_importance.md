# Compute Permutation Feature Importance

Calculates permutation-based feature importance using DALEX.

## Usage

``` r
xai_importance(
  explainer,
  n_perm = 10,
  loss_function = "1 - auc",
  features = NULL
)
```

## Arguments

- explainer:

  A DALEX explainer

- n_perm:

  Number of permutations (default: 10)

- loss_function:

  Loss function ("1-auc" for classification)

- features:

  Optional: subset of features to evaluate

## Value

A model_parts object with importance scores
