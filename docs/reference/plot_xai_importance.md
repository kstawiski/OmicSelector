# Plot XAI Feature Importance

Creates a feature importance plot with correlation warnings.

## Usage

``` r
plot_xai_importance(xai_results, top_n = 15, show_cor_warning = TRUE)
```

## Arguments

- xai_results:

  Output from xai_pipeline

- top_n:

  Number of top features to show (default: 15)

- show_cor_warning:

  Highlight correlated features (default: TRUE)

## Value

A ggplot object
