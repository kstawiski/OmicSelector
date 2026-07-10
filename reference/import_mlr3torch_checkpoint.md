# Import mlr3torch Checkpoint

Loads a torch checkpoint into a trained learner for inference, or
attaches it for potential fine-tuning if the learner is not yet trained.

## Usage

``` r
import_mlr3torch_checkpoint(learner, path, strict = TRUE)
```

## Arguments

- learner:

  mlr3 learner or GraphLearner.

- path:

  Checkpoint path created by
  [`export_mlr3torch_checkpoint()`](https://kstawiski.github.io/OmicSelector/reference/export_mlr3torch_checkpoint.md).

- strict:

  Logical, enforce exact layer matching.

## Value

The updated learner.
