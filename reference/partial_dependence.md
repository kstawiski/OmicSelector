# Compute Partial Dependence

Computes partial dependence profiles showing how a feature affects
predictions, marginalizing over other features.

## Usage

``` r
partial_dependence(explainer, features, n_grid = 20)
```

## Arguments

- explainer:

  An OmicExplainer object

- features:

  Character vector of feature names to analyze

- n_grid:

  Number of grid points (default: 20)

## Value

A data frame with partial dependence values
