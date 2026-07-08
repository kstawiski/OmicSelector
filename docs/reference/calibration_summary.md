# Calibration Summary for Model Results

Computes comprehensive calibration metrics from benchmark results.
Intended for use with BenchmarkService outputs.

This function is metrics-only and does NOT fit calibration repair
models. To fit calibrators, use the separate functions
[`fit_platt_scaling()`](https://kstawiski.github.io/OmicSelector/reference/fit_platt_scaling.md),
[`fit_isotonic_calibration()`](https://kstawiski.github.io/OmicSelector/reference/fit_isotonic_calibration.md),
or
[`fit_temperature_scaling()`](https://kstawiski.github.io/OmicSelector/reference/fit_temperature_scaling.md)
on a held-out calibration set (NOT the evaluation/test set).

## Usage

``` r
calibration_summary(probs, labels)
```

## Arguments

- probs:

  Vector or list of predicted probabilities

- labels:

  Vector or list of true labels

## Value

A CalibrationResult object with calibration metrics

## Details

To avoid data leakage, calibration repair models should be fit on a
separate calibration set, not on the evaluation data used for metrics.
This function only computes metrics; fitting calibrators is a separate
step.

## See also

[`fit_platt_scaling`](https://kstawiski.github.io/OmicSelector/reference/fit_platt_scaling.md),
[`fit_isotonic_calibration`](https://kstawiski.github.io/OmicSelector/reference/fit_isotonic_calibration.md),
[`fit_temperature_scaling`](https://kstawiski.github.io/OmicSelector/reference/fit_temperature_scaling.md)
