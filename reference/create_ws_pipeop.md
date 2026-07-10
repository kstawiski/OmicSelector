# Apply Within-Sample Normalization to OmicPipeline

Creates a PipeOp for within-sample normalization that can be integrated
into the OmicSelector mlr3 pipeline. This ensures within-sample
normalization is applied inside CV folds with zero leakage.

## Usage

``` r
create_ws_pipeop(method = "minmax", id = "ws_norm")
```

## Arguments

- method:

  One of "minmax", "rank", "zscore", "logratio", "ratio_image"

- id:

  PipeOp identifier

## Value

A PipeOpTaskPreproc for within-sample normalization
