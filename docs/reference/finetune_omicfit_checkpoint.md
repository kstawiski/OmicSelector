# Fine-tune OmicFit from mlr3torch Checkpoint

Convenience wrapper that fine-tunes an OmicFit learner and returns a new
OmicFit object.

## Usage

``` r
finetune_omicfit_checkpoint(fit, task, path, epochs = NULL,
  seed = NULL, device = NULL, strict = TRUE)
```

## Arguments

- fit:

  OmicFit object returned by `OmicPipeline$fit()`.

- task:

  mlr3 Task to train on.

- path:

  Checkpoint path.

- epochs:

  Optional number of epochs to set before training.

- seed:

  Optional seed.

- device:

  Optional device ("cpu" or "cuda") if supported by learner.

- strict:

  Logical, enforce exact layer matching.

## Value

A new OmicFit object.
