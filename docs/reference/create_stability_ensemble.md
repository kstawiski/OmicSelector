# Create a Stability Ensemble with Presets

Convenience function to create a StabilityEnsemble with common presets.

## Usage

``` r
create_stability_ensemble(preset = "default", ...)
```

## Arguments

- preset:

  One of "default", "fast", "thorough"

- ...:

  Additional arguments passed to StabilityEnsemble\$new()

## Value

A StabilityEnsemble object

## Examples

``` r
if (FALSE) { # \dontrun{
# Default: 100 bootstraps, glmnet base learner
ensemble <- create_stability_ensemble("default")

# Fast: 50 bootstraps for quick exploration
ensemble <- create_stability_ensemble("fast")

# Thorough: 500 bootstraps for publication-quality results
ensemble <- create_stability_ensemble("thorough")

# Fit and predict
ensemble$fit(task)
probs <- ensemble$predict(task)
} # }
```
