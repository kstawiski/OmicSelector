# Score a Direction-Split Panel by Within-Sample Ranks

Computes a singscore-style score from within-sample fractional ranks:
mean rank of up features minus mean rank of down features. The split
cardinality is carried in the return attributes for null-QC checks.

## Usage

``` r
os_singscore(X, up, down)
```

## Arguments

- X:

  Numeric matrix or data frame.

- up:

  Character vector of up-in-positive feature names.

- down:

  Character vector of down-in-positive feature names.

## Value

Numeric score vector.
