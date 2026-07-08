# Aggregate SHAP Values from List

Aggregates SHAP results into a matrix.

## Usage

``` r
aggregate_shap_values(shap_results, feature_names)
```

## Arguments

- shap_results:

  List of SHAP data.frames

- feature_names:

  Character vector of feature names

## Value

Matrix of SHAP values (samples x features)
