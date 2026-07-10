# Mahalanobis Anomaly Score from a Reference Set

Mahalanobis Anomaly Score from a Reference Set

## Usage

``` r
os_mahalanobis_score(X, reference_rows, ridge = 1e-06)
```

## Arguments

- X:

  Numeric matrix or data frame.

- reference_rows:

  Integer, logical, or character row selector for the reference set.

- ridge:

  Non-negative diagonal ridge added to the covariance matrix.

## Value

Numeric squared Mahalanobis distance per row.
