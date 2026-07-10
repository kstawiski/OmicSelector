# Fine-tune from mlr3torch Checkpoint

Loads a checkpoint and trains the learner on a new task, attempting to
warm-start when supported by the mlr3torch learner.

## Usage

``` r
finetune_mlr3torch_checkpoint(learner, task, path, epochs = NULL,
  seed = NULL, device = NULL, strict = TRUE)
```

## Arguments

- learner:

  mlr3 learner or GraphLearner.

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

Trained learner.
