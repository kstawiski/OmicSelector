# Make MLP Learner (Alias)

Convenient alias for create_mlp_learner() matching the naming convention
of other factory functions (make_autotuner\_\*).

## Usage

``` r
make_mlp_learner(
  hidden_layers = c(64, 32),
  dropout = 0.2,
  activation = c("relu", "tanh", "sigmoid"),
  epochs = 100,
  batch_size = 32,
  lr = 0.001,
  early_stopping = TRUE,
  patience = 10,
  seed = NULL,
  device = "cpu"
)
```

## Arguments

- hidden_layers:

  Integer vector specifying neurons per hidden layer. Default: c(64, 32)
  for a 2-layer network.

- dropout:

  Dropout rate between 0 and 1. Default: 0.2.

- activation:

  Activation function: "relu" (default), "tanh", "sigmoid".

- epochs:

  Maximum training epochs. Default: 100.

- batch_size:

  Training batch size. Default: 32.

- lr:

  Learning rate. Default: 0.001.

- early_stopping:

  Logical, use early stopping. Default: TRUE.

- patience:

  Early stopping patience (epochs). Default: 10.

- seed:

  Random seed for reproducibility.

- device:

  "cpu" or "cuda" for GPU. Default: "cpu".

## Value

An mlr3 Learner object
