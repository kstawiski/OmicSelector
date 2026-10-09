# Platt Scaling (Logistic Calibration)

Fits logistic regression on the clipped probability score to recalibrate
probabilities. The input is a probability in \[0, 1\], not a log-odds
value. Parameters should be estimated on a held-out calibration set.

## Usage

``` r
fit_platt_scaling(probs, labels)
```

## Arguments

- probs:

  Calibration set predicted probabilities in \[0, 1\]

- labels:

  Calibration set labels

## Value

A function that recalibrates probabilities

## Details

Fits `glm(y ~ score, family = binomial())`, where `score` is the
supplied probability clipped to `[1e-10, 1 - 1e-10]`. Predictions use
`plogis(intercept + slope * clipped_score)` with the same clipping. No
log-odds transformation is applied to the input.

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
