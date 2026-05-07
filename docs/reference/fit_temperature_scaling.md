# Temperature Scaling

Fits temperature scaling for probability calibration. Uses a single
temperature parameter T to soften/sharpen predictions.

## Usage

``` r
fit_temperature_scaling(probs, labels)
```

## Arguments

- probs:

  Calibration set predicted probabilities

- labels:

  Calibration set labels

## Value

A function that recalibrates probabilities

## Details

Temperature scaling modifies the logits by dividing by T:
\$\$p\_{calibrated} = \sigma(logit(p) / T)\$\$

\- T \> 1: softens predictions (moves toward 0.5) - T \< 1: sharpens
predictions (moves toward 0 or 1) - T = 1: no change

Commonly used for neural network calibration.
