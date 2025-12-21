# Compute Permutation Feature Importance

Computes model-agnostic feature importance by permuting features and
measuring the drop in performance.

## Usage

``` r
feature_importance(
  explainer,
  loss_function = c("auc_loss", "accuracy_loss", "brier"),
  n_permutations = 10,
  parallel = FALSE
)
```

## Arguments

- explainer:

  An OmicExplainer object

- loss_function:

  Loss function to use (default: "1 - AUC")

- n_permutations:

  Number of permutations per feature (default: 10)

- parallel:

  Logical, use parallel computation (default: FALSE)

## Value

A data frame with feature importance scores

## Details

Permutation importance measures how much model performance degrades when
a feature's values are randomly shuffled. Higher values indicate more
important features.

For correlated features, importance may be shared or masked. Use
\`feature_importance_groups()\` for correlated feature handling.

## Examples

``` r
if (FALSE) { # \dontrun{
importance <- feature_importance(explainer, n_permutations = 20)
head(importance[order(-importance$importance), ])
} # }
```
