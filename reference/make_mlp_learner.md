# Create Omics-Optimized MLP Learner

Creates a Multi-Layer Perceptron with regularization optimized for
high-dimensional omics data.

## Usage

``` r
make_mlp_learner(
  n_hidden = 128L,
  n_layers = 2L,
  dropout = 0.5,
  batch_size = 32L,
  epochs = 100L,
  learning_rate = 0.001,
  weight_decay = 0.01,
  early_stopping_patience = 10L
)
```

## Arguments

- n_hidden:

  Number of hidden units per layer (default: 128)

- n_layers:

  Number of hidden layers (default: 2)

- dropout:

  Dropout rate (default: 0.5)

- batch_size:

  Training batch size (default: 32)

- epochs:

  Number of training epochs (default: 100)

- learning_rate:

  Learning rate (default: 0.001)

- weight_decay:

  L2 regularization (default: 0.01)

- early_stopping_patience:

  Patience for early stopping (default: 10)

## Value

An mlr3 Learner object

## Details

The MLP is configured with: - High dropout (0.5) to prevent overfitting
on high-dimensional data - Strong L2 regularization - Batch
normalization between layers - Adam optimizer with early stopping

## Examples

``` r
if (FALSE) { # \dontrun{
learner <- make_mlp_learner(n_hidden = 64, dropout = 0.5)
learner$train(task)
} # }
```
