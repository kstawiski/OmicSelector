# Import mlr3torch Checkpoint into OmicFit

Convenience wrapper that loads a checkpoint into the learner inside an
OmicFit object.

## Usage

``` r
import_omicfit_checkpoint(fit, path, strict = TRUE)
```

## Arguments

- fit:

  OmicFit object returned by `OmicPipeline$fit()`.

- path:

  Checkpoint path created by
  [`export_mlr3torch_checkpoint()`](https://kstawiski.github.io/OmicSelector/reference/export_mlr3torch_checkpoint.md).

- strict:

  Logical, enforce exact layer matching.

## Value

Updated OmicFit object.
