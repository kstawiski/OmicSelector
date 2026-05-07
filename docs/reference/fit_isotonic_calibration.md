# Isotonic Regression Calibration

Fits isotonic (monotonic) regression for probability calibration.
Non-parametric method that produces a monotonically increasing
calibration mapping.

## Usage

``` r
fit_isotonic_calibration(probs, labels)
```

## Arguments

- probs:

  Calibration set predicted probabilities

- labels:

  Calibration set labels

## Value

A function that recalibrates probabilities

## Details

Isotonic regression finds the best monotonically increasing function
that maps scores to calibrated probabilities. It's more flexible than
Platt scaling but requires more calibration data.

IMPORTANT: Must be fit on a held-out calibration set, not training data,
to avoid leakage.

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit on calibration set
calibrator <- fit_isotonic_calibration(cal_probs, cal_labels)

# Apply to new predictions
calibrated <- calibrator(new_probs)
} # }
```
