# Compute Apparent or Cross-Fitted Calibrated Brier Score

Compute Apparent or Cross-Fitted Calibrated Brier Score

## Usage

``` r
os_calibrated_brier(y, score, folds = NULL, positive = NULL)
```

## Arguments

- y:

  Binary outcome vector.

- score:

  Numeric prediction score.

- folds:

  Optional list of held-out row indices for cross-fitting Platt
  calibration.

- positive:

  Optional positive-class label.

## Value

An `os_calibrated_brier` object.
