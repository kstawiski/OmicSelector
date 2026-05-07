# Run Complete XAI Pipeline

Runs the full explainability pipeline: importance, PDPs, SHAP,
correlations.

## Usage

``` r
xai_pipeline(
  learner,
  task,
  top_k = 10,
  n_shap_obs = 5,
  cor_threshold = 0.7,
  features = NULL,
  verbose = TRUE
)
```

## Arguments

- learner:

  A trained mlr3 learner

- task:

  The mlr3 task

- top_k:

  Number of top features to analyze in detail (default: 10)

- n_shap_obs:

  Number of observations for SHAP (default: 5)

- cor_threshold:

  Correlation threshold for warnings (default: 0.7)

- features:

  Optional: pre-selected feature subset

- verbose:

  Print progress (default: TRUE)

## Value

A list containing all XAI results

## Examples

``` r
if (FALSE) { # \dontrun{
results <- xai_pipeline(learner, task, top_k = 10)
plot(results$importance)
print(results$correlations$high_cor_pairs)
} # }
```
