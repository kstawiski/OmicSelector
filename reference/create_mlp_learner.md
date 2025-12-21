# Create MLP Learner via mlr3torch

Creates a Multi-Layer Perceptron (MLP) learner for classification tasks.
Uses mlr3torch for native torch integration within the mlr3 ecosystem.

## Usage

``` r
create_mlp_learner(
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

An mlr3 Learner object (or a fallback if mlr3torch unavailable)

## Details

The MLP architecture: - Input layer: matches number of features - Hidden
layers: as specified, with batch normalization and dropout - Output
layer: 2 neurons (binary classification) with softmax

If mlr3torch is not installed, returns a fallback learner (glmnet) with
a warning.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic MLP
mlp <- create_mlp_learner()

# Deeper network with more regularization
mlp_deep <- create_mlp_learner(
  hidden_layers = c(128, 64, 32),
  dropout = 0.3,
  epochs = 200
)

# Use in OmicPipeline
pipeline <- OmicPipeline$new(data, target = "outcome")
learner <- pipeline$create_graph_learner(
  filter = "anova",
  model = mlp,  # Pass the MLP learner
  n_features = 50
)
} # }
```
