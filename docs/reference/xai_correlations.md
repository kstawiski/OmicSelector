# Compute Correlation Diagnostics for Features

Identifies highly correlated feature pairs that may cause misleading
SHAP values or permutation importance.

## Usage

``` r
xai_correlations(data, features = NULL, method = "spearman", threshold = 0.7)
```

## Arguments

- data:

  Data frame of features

- features:

  Optional: subset of features to analyze

- method:

  Correlation method ("spearman", "pearson")

- threshold:

  Correlation threshold for warnings (default: 0.7)

## Value

A list with correlation matrix, clusters, and warnings
