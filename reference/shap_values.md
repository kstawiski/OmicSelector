# Compute SHAP-like Values

Computes Shapley-like additive feature attributions for individual
predictions. Uses kernel SHAP approximation for efficiency.

## Usage

``` r
shap_values(explainer, new_data, n_samples = 100)
```

## Arguments

- explainer:

  An OmicExplainer object

- new_data:

  Data for which to compute explanations

- n_samples:

  Number of samples for approximation (default: 100)

## Value

A matrix of SHAP values (instances x features)

## Details

SHAP (SHapley Additive exPlanations) values explain individual
predictions as the sum of feature contributions. For each prediction:

prediction = base_value + sum(SHAP_values)

Positive SHAP value = feature pushes prediction higher Negative SHAP
value = feature pushes prediction lower
