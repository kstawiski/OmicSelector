# Create an FT-Transformer Learner

Create an \`mlr3\` classification learner backed by
\`rtdl_revisiting_models::FTTransformer\` through \`reticulate\`.

## Usage

``` r
make_fttransformer_learner(
  n_blocks = 3L,
  d_block = 192L,
  attention_dropout = 0.2,
  ffn_dropout = 0.2,
  ffn_d_hidden = NULL,
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

make_tabtransformer_learner(...)
```

## Arguments

- n_blocks:

  Number of transformer blocks. Default: \`3\`.

- d_block:

  Transformer width. Default: \`192\`.

- attention_dropout:

  Attention dropout. Default: \`0.20\`.

- ffn_dropout:

  Feed-forward dropout. Default: \`0.20\`.

- ffn_d_hidden:

  Hidden width in the feed-forward network. Default: \`NULL\`, which
  leaves the backend default unchanged.

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

- ...:

  Additional arguments forwarded to \[make_fttransformer_learner()\].

## Value

An \`mlr3::LearnerClassif\` with \`predict_type = "prob"\`.

## References

Gorishniy, Y., et al. (2021). Revisiting Deep Learning Models for
Tabular Data.
