# Export Model to ONNX Format

Exports the estimator component of an mlr3 learner to ONNX format for
cross-language deployment. Supports xgboost and glmnet learners.

## Usage

``` r
export_onnx(learner, path, task, export_preprocessing = TRUE)
```

## Arguments

- learner:

  A trained mlr3 Learner or GraphLearner

- path:

  File path for the ONNX model (with .onnx extension)

- task:

  The mlr3 Task used for training (for metadata)

- export_preprocessing:

  Logical, whether to also export preprocessing manifest

## Value

List with export status and paths to exported files

## Details

ONNX export is currently supported for: - xgboost: via native xgboost
save + Python conversion (most reliable) - glmnet: via coefficient
extraction (linear models)

For xgboost, this function exports the native .xgb format. A separate
Python script can convert to ONNX. For glmnet, coefficients are exported
as JSON.

For other learners, only the preprocessing manifest is exported. Use
vetiver for R-based deployment of these models.

## Examples

``` r
if (FALSE) { # \dontrun{
# Train xgboost model
task <- tsk("iris")
learner <- lrn("classif.xgboost", nrounds = 50)
learner$train(task)

# Export to ONNX-compatible format
export_onnx(learner, "model.onnx", task)
} # }
```
