# Compute SHAP Values for Observations

Computes SHAP values for individual predictions.

## Usage

``` r
xai_shap(explainer, new_observations, B = 25)
```

## Arguments

- explainer:

  A DALEX explainer

- new_observations:

  Data frame of observations to explain

- B:

  Number of random orderings (default: 25)

## Value

A list of predict_parts objects (one per observation)
