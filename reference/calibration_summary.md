# Calibration Summary for Model Results

Computes comprehensive calibration metrics from benchmark results.
Intended for use with BenchmarkService outputs.

## Usage

``` r
calibration_summary(probs, labels, repair = FALSE)
```

## Arguments

- probs:

  Vector or list of predicted probabilities

- labels:

  Vector or list of true labels

- repair:

  Logical, whether to fit calibration repair methods

## Value

A CalibrationResult object with metrics and optional repairers
