# Get Reliable SHAP Features

Returns only the features whose SHAP values are reliable (no high
correlations).

## Usage

``` r
get_reliable_shap_features(shap_result, top_n = NULL)
```

## Arguments

- shap_result:

  Result from compute_shap_with_warnings()

- top_n:

  Number of top features to return (default: all)

## Value

Named vector of feature importance for reliable features only
