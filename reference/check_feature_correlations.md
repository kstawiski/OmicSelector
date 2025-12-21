# Check Feature Correlations

Identifies highly correlated features that may lead to misleading
importance scores. When features are correlated, importance can be
shared, masked, or unstable.

## Usage

``` r
check_feature_correlations(
  data,
  threshold = 0.7,
  method = c("spearman", "pearson")
)
```

## Arguments

- data:

  Data frame with features

- threshold:

  Correlation threshold (default: 0.7)

- method:

  Correlation method: "pearson", "spearman" (default: "spearman")

## Value

A list with correlation warnings

## Details

High correlation between features causes problems for interpretation: -
Permutation importance is shared between correlated features - SHAP
values become unstable - One feature can "mask" another

This function identifies correlated pairs and suggests grouped analysis.
