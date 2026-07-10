# Create a Hybrid Sequential Feature Selector

Convenience function to create a SequentialSelector with common presets.

## Usage

``` r
create_hsfs_selector(preset = "default", ...)
```

## Arguments

- preset:

  One of "default", "aggressive", "conservative", or "custom"

- ...:

  Additional arguments passed to SequentialSelector\$new()

## Value

A SequentialSelector object

## Examples

``` r
if (FALSE) { # \dontrun{
# Default pipeline: variance -> ANOVA(5000) -> RF_importance(1000) -> AUC_filter
selector <- create_hsfs_selector("default")

# Aggressive reduction for very high-dimensional data
selector <- create_hsfs_selector("aggressive")

# Conservative - keep more features
selector <- create_hsfs_selector("conservative")

# Create GraphLearner
learner <- selector$create_learner(model = "xgboost")
} # }
```
