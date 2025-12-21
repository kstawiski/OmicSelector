# Create Complete Export Bundle

Creates a complete export bundle containing: - Vetiver model for R
deployment - ONNX-compatible model files (if supported) - Preprocessing
manifest - Feature metadata

## Usage

``` r
export_bundle(
  learner,
  model_name,
  task,
  output_dir,
  description = NULL,
  create_api = FALSE
)
```

## Arguments

- learner:

  A trained mlr3 Learner or GraphLearner

- model_name:

  Name for the model bundle

- task:

  The mlr3 Task used for training

- output_dir:

  Directory for export files

- description:

  Optional description

- create_api:

  Logical, whether to generate Plumber API script

## Value

List with paths to all exported files

## Examples

``` r
if (FALSE) { # \dontrun{
# Create complete export bundle
bundle <- export_bundle(
  learner = trained_learner,
  model_name = "biomarker_classifier",
  task = my_task,
  output_dir = "model_export",
  create_api = TRUE
)
} # }
```
