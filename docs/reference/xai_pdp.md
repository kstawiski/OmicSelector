# Compute Partial Dependence Plots

Computes PDPs for specified features using DALEX.

## Usage

``` r
xai_pdp(explainer, features, grid_n = 20, grid_type = "quantile")
```

## Arguments

- explainer:

  A DALEX explainer

- features:

  Features to compute PDPs for

- grid_n:

  Number of grid points (default: 20)

- grid_type:

  "quantile" or "uniform"

## Value

A model_profile object
