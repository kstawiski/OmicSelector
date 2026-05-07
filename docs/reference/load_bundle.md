# Load Exported Model Bundle

Loads a model from an export bundle created by export_bundle().

## Usage

``` r
load_bundle(bundle_dir, model_name = NULL)
```

## Arguments

- bundle_dir:

  Directory containing the export bundle

- model_name:

  Name of the model (defaults to directory name)

## Value

List with learner and task_info
