# Validate Synthetic Data Quality

Computes quality metrics for synthetic data compared to real data.

## Usage

``` r
validate_synthetic(real, synthetic, target, features = NULL)
```

## Arguments

- real:

  Real data (data.frame)

- synthetic:

  Synthetic data (data.frame)

- target:

  Target column name

- features:

  Feature columns to evaluate (default: all numeric)

## Value

A list with quality metrics

## Details

Metrics computed: - \*\*Distribution similarity\*\*: KS statistic per
feature - \*\*Correlation preservation\*\*: Correlation matrix
similarity - \*\*Class balance\*\*: Proportion comparison - \*\*Outlier
rate\*\*: Proportion of extreme values
