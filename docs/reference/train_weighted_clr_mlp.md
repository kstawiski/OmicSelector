# Train a Weighted CLR + MLP Classifier

Fits a hemolysis-aware weighted-CLR analogue of \[train_clr_mlp()\]
using the same backends and normalization scheme as the dense CoDA
baseline.

## Usage

``` r
train_weighted_clr_mlp(
  X_train,
  y_train,
  weights = NULL,
  sensitivity = .paper1_hemolysis_coefficients(),
  shrink_features = TRUE,
  epochs = 30L,
  lr = 0.001,
  batch_size = 256L,
  class_weight = TRUE,
  device = NULL,
  verbose = TRUE,
  seed = 42L,
  positive = NULL,
  input_scale = c("log", "raw"),
  backend = c("auto", "torch", "nnet", "glm"),
  hidden_units = c(152L, 136L),
  dropout = 0.3
)
```

## Arguments

- X_train:

  Numeric training matrix with samples in rows and features in columns.

- y_train:

  Binary outcome vector.

- weights:

  Optional positive numeric vector of feature weights.

- sensitivity:

  Optional named positive numeric vector of hemolysis coefficients used
  when \`weights = NULL\`.

- shrink_features:

  Logical; whether to shrink centered coordinates by feature weight
  after weighted centering.

- epochs:

  Number of training epochs. Defaults to \`30\`.

- lr:

  Learning rate. Defaults to \`0.001\`.

- batch_size:

  Mini-batch size. Defaults to \`256\`.

- class_weight:

  Logical; if \`TRUE\`, use inverse-frequency weighting.

- device:

  One of \`"cpu"\`, \`"cuda"\`, or \`NULL\`/\`"auto"\`.

- verbose:

  Logical; print training progress. Defaults to \`TRUE\`.

- seed:

  Integer random seed. Defaults to \`42\`.

- positive:

  Optional positive-class label for non-numeric outcomes.

- input_scale:

  Either \`"log"\` (default) or \`"raw"\`.

- backend:

  One of \`"auto"\`, \`"torch"\`, \`"nnet"\`, or \`"glm"\`.

- hidden_units:

  Integer vector with the two hidden-layer sizes used by the torch
  backend. Defaults to \`c(152, 136)\`.

- dropout:

  Dropout rate used by the torch backend. Defaults to \`0.3\`.

## Value

An object of class \`"omicselector_weighted_clr_mlp"\`.
