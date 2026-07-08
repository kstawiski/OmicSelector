# Within-Sample Min-Max Normalization

Scales each sample's feature values to \[0,1\] range using only that
sample's minimum and maximum. Completely batch-immune.

## Usage

``` r
ws_minmax(x)
```

## Arguments

- x:

  A numeric matrix (samples x features) or vector (single sample)

## Value

Matrix of same dimensions with values in \[0,1\]

## Examples

``` r
if (FALSE) { # \dontrun{
# Single sample
ws_minmax(c(5.2, 3.1, 8.7, 2.4))
# Matrix (samples x features)
ws_minmax(matrix(rnorm(100), nrow=10))
} # }
```
