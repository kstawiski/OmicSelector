# Compute conformal anomaly p-value for new samples

Scores one or more new samples against the calibrated conformal anomaly
model from
[`fit_conformal_anomaly`](https://kstawiski.github.io/OmicSelector/reference/fit_conformal_anomaly.md).
The conformal p-value is the fraction of calibration scores at least as
extreme as the test score: \$\$p = \frac{1 + \\\\s\_{\mathrm{cal}} \geq
s\_{\mathrm{test}}\\}{1 + n\_{\mathrm{cal}}}.\$\$

## Usage

``` r
os_conformal_anomaly(x, fit, alpha = 0.05)
```

## Arguments

- x:

  Numeric vector or matrix (samples \\\times\\ features).

- fit:

  A `conformal_anomaly_fit` object from
  [`fit_conformal_anomaly`](https://kstawiski.github.io/OmicSelector/reference/fit_conformal_anomaly.md).

- alpha:

  Significance level for the binary anomaly call. Default 0.05.

## Value

List with `score` (raw k-NN distance), `p_value` (conformal p in
\\\[0,1\]\\), `is_anomaly` (logical at level `alpha`), `alpha`, and
`calibration_summary`.
