# Within-Sample Rank Normalization

Replaces each sample's feature values with their within-sample
fractional ranks. Completely ordinal and batch-immune.

## Usage

``` r
ws_rank(x)
```

## Arguments

- x:

  A numeric matrix (samples x features) or vector

## Value

Matrix with values in (0,1\] representing fractional ranks
