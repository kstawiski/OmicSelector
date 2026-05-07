# Compute Binary Clinical Operating Points

Summarizes sensitivity and specificity at prespecified targets. Larger
scores are treated as more positive.

## Usage

``` r
os_operating_points(
  y,
  score,
  specificity_targets = c(0.9, 0.95, 0.98),
  sensitivity_targets = c(0.9, 0.95),
  positive = NULL
)
```

## Arguments

- y:

  Binary outcome vector. The second factor level is treated as positive
  unless `positive` is supplied.

- score:

  Numeric prediction score.

- specificity_targets:

  Specificity thresholds to evaluate.

- sensitivity_targets:

  Sensitivity thresholds to evaluate.

- positive:

  Optional positive-class label.

## Value

An `os_operating_points` object.
