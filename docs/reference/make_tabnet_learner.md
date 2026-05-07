# Create a TabNet Learner

Create an \`mlr3\` classification learner backed by
\`pytorch_tabnet::TabNetClassifier\` through \`reticulate\`.

## Usage

``` r
make_tabnet_learner(
  n_d = 16L,
  n_a = 16L,
  n_steps = 4L,
  gamma = 1.5,
  n_independent = 2L,
  n_shared = 2L,
  lambda_sparse = 1e-04,
  batch_size = 128L,
  virtual_batch_size = 32L,
  max_epochs = 200L,
  patience = 20L,
  learning_rate = 0.02,
  weight_decay = 1e-05,
  mask_type = c("entmax", "sparsemax"),
  validation_fraction = 0.15,
  class_weights = "balanced",
  device = c("auto", "cpu", "cuda"),
  seed = 20260331L
)
```

## Arguments

- n_d:

  Width of the decision layer. Default: \`16\`.

- n_a:

  Width of the attention layer. Default: \`16\`.

- n_steps:

  Number of sequential attention steps. Default: \`4\`.

- gamma:

  Feature reusage penalty. Default: \`1.5\`.

- n_independent:

  Number of independent GLU blocks. Default: \`2\`.

- n_shared:

  Number of shared GLU blocks. Default: \`2\`.

- lambda_sparse:

  Sparsity regularization strength. Default: \`1e-4\`.

- batch_size:

  Minibatch size. Default: \`128\`.

- virtual_batch_size:

  Ghost batch normalization size. Default: \`32\`.

- max_epochs:

  Maximum epochs. Default: \`200\`.

- patience:

  Early stopping patience. Default: \`20\`.

- learning_rate:

  Optimizer learning rate. Default: \`0.02\`.

- weight_decay:

  L2 regularization for Adam. Default: \`1e-5\`.

- mask_type:

  Mask function, \`"sparsemax"\` or \`"entmax"\`. Default: \`"entmax"\`.

- validation_fraction:

  Fraction of the training fold used for internal early stopping.
  Default: \`0.15\`.

- class_weights:

  Class weighting scheme. Use \`"balanced"\` (default), \`"none"\`, or a
  numeric vector with one weight per class. Named vectors are matched
  against the task class labels.

- device:

  \`"auto"\` (default), \`"cpu"\`, or \`"cuda"\`.

- seed:

  Random seed used inside the backend. Default: \`20260331\`.

## Value

An \`mlr3::LearnerClassif\` with \`predict_type = "prob"\`.

## References

Arik, S. O., & Pfister, T. (2021). TabNet: Attentive Interpretable
Tabular Learning.
