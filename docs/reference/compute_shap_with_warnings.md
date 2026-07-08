# Compute SHAP Values with Correlation Warnings

Computes SHAP values for model predictions and flags features that have
high correlations with other features, which can make SHAP
interpretations unreliable.

## Usage

``` r
compute_shap_with_warnings(
  model,
  data,
  n_samples = 100,
  correlation_threshold = 0.7,
  method = "permutation",
  verbose = TRUE
)
```

## Arguments

- model:

  A trained model (mlr3 Learner, caret model, or predict function)

- data:

  Data.frame of features used for SHAP computation

- n_samples:

  Number of samples for SHAP approximation (default: 100)

- correlation_threshold:

  Correlation threshold for warnings (default: 0.7)

- method:

  SHAP method: "permutation" or "kernel" (default: "permutation")

- verbose:

  Print progress messages

## Value

A list with: - shap_values: Matrix of SHAP values (samples x features) -
feature_importance: Mean absolute SHAP per feature - correlations:
Correlation matrix of features - warnings: Data.frame of correlated
feature pairs with warnings - reliable_features: Features with no high
correlations (safe to interpret)

## Examples

``` r
if (FALSE) { # \dontrun{
# Compute SHAP with correlation warnings
result <- compute_shap_with_warnings(model, test_data)

# Check which features have reliable SHAP interpretations
print(result$reliable_features)

# View correlation warnings
print(result$warnings)
} # }
```
