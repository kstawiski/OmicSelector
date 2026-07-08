# Platt Scaling (Logistic Calibration)

Fits a logistic regression model to recalibrate probabilities.
Parameters should be estimated on a held-out calibration set.

## Usage

``` r
fit_platt_scaling(probs, labels)
```

## Arguments

- probs:

  Calibration set predicted probabilities

- labels:

  Calibration set labels

## Value

A function that recalibrates probabilities

## Details

Platt scaling fits: P(y=1\|f) = 1 / (1 + exp(A\*f + B)) where f is the
original score/probability and A, B are fitted parameters.

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit on calibration set
cal_probs <- runif(100)
cal_labels <- rbinom(100, 1, cal_probs)
calibrator <- fit_platt_scaling(cal_probs, cal_labels)

# Apply to new predictions
new_probs <- c(0.2, 0.5, 0.8)
calibrated <- calibrator(new_probs)
} # }
```
