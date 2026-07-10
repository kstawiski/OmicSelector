# Export mlr3torch Checkpoint from OmicFit

Convenience wrapper that exports a checkpoint from an OmicFit object.

## Usage

``` r
export_omicfit_checkpoint(fit, path, include_optimizer = FALSE,
  metadata = NULL)
```

## Arguments

- fit:

  OmicFit object returned by `OmicPipeline$fit()`.

- path:

  Output file path (recommended: .pt).

- include_optimizer:

  Logical; attempts to save optimizer state if available.

- metadata:

  Optional list of metadata to store alongside weights.

## Value

The file path (invisibly).
