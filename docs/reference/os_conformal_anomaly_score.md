# Conformal Healthy-Reference Anomaly Score

Computes a split-conformal anomaly score using Mahalanobis distances
from reference-training controls and calibration controls. The returned
score is `1 - conformal p-value`; larger values are more anomalous.

## Usage

``` r
os_conformal_anomaly_score(X, train_rows, calibration_rows, ridge = 1e-06)
```

## Arguments

- X:

  Numeric matrix or data frame.

- train_rows:

  Reference-training row selector.

- calibration_rows:

  Calibration row selector.

- ridge:

  Non-negative diagonal ridge for covariance.

## Value

Numeric score vector.
