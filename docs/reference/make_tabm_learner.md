# Create a TabM Learner

Create an \`mlr3\` classification learner for the modern MLP baseline
TabM through \`reticulate\`. The wrapper accepts either the official
\`tabm\` package or a compatible \`rtdl_revisiting_models\` build that
exposes \`TabM\`.

## Usage

``` r
make_tabm_learner(
  d_block = 128L,
  n_blocks = 3L,
  dropout = 0.1,
  k = NULL,
  batch_size = 128L,
  max_epochs = 200L,
  patience = 20L,
  learning_rate = 0.001,
  weight_decay = 1e-05,
  validation_fraction = 0.15,
  class_weights = "balanced",
  device = c("auto", "cpu", "cuda"),
  seed = 20260331L
)
```

## Arguments

- d_block:

  Hidden block width. Default: \`128\`.

- n_blocks:

  Number of MLP blocks. Default: \`3\`.

- dropout:

  Model dropout. Default: \`0.10\`.

- k:

  Optional TabM ensemble width if supported by the backend. Default:
  \`NULL\`.

- batch_size:

  Minibatch size. Default: \`128\`.

- max_epochs:

  Maximum epochs. Default: \`200\`.

- patience:

  Early stopping patience. Default: \`20\`.

- learning_rate:

  Optimizer learning rate. Default: \`1e-3\`.

- weight_decay:

  AdamW weight decay. Default: \`1e-5\`.

- validation_fraction:

  Fraction of the training fold used for internal validation. Default:
  \`0.15\`.

- class_weights:

  Class weighting scheme. Use \`"balanced"\` (default), \`"none"\`, or a
  numeric vector with one weight per class.

- device:

  \`"auto"\` (default), \`"cpu"\`, or \`"cuda"\`.

- seed:

  Random seed used inside the backend. Default: \`20260331\`.

## Value

An \`mlr3::LearnerClassif\` with \`predict_type = "prob"\`.

## References

Gorishniy, Y., et al. (2024). TabM: Advancing Tabular Deep Learning with
Parameter-Efficient Ensembling.
