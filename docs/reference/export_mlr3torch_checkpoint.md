# Export mlr3torch Checkpoint

Saves a trained mlr3torch learner's torch module state to disk.

## Usage

``` r
export_mlr3torch_checkpoint(learner, path, include_optimizer = FALSE,
  metadata = NULL)
```

## Arguments

- learner:

  Trained mlr3 learner or GraphLearner.

- path:

  Output file path (recommended: .pt).

- include_optimizer:

  Logical; attempts to save optimizer state if available.

- metadata:

  Optional list of metadata to store alongside weights.

## Value

The file path (invisibly).
