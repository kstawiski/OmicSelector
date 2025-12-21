# Generate Correlation Warnings

Creates warning messages for correlated feature pairs.

## Usage

``` r
generate_correlation_warnings(high_cor_pairs, feature_importance, threshold)
```

## Arguments

- high_cor_pairs:

  Data.frame of correlated pairs

- feature_importance:

  Named vector of feature importance

- threshold:

  Correlation threshold

## Value

Data.frame with warnings
