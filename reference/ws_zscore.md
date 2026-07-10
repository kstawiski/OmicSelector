# Within-Sample Z-Score Normalization

Standardizes each sample across its own features (mean=0, sd=1 within
sample). Batch-immune because it uses only within-sample statistics.

## Usage

``` r
ws_zscore(x)
```

## Arguments

- x:

  A numeric matrix (samples x features) or vector

## Value

Matrix with within-sample standardized values
