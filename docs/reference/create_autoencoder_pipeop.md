# Create Autoencoder PipeOp

Creates an mlr3pipelines PipeOp that trains a torch autoencoder on the
training data and replaces features with latent representations.

## Usage

``` r
create_autoencoder_pipeop(latent_dim = 32L, hidden_layers = c(128L, 64L),
  activation = "relu", dropout = 0.1, epochs = 100L, batch_size = 32L,
  lr = 0.001, weight_decay = 0, loss = "mse", validation_split = 0.1,
  early_stopping = TRUE, patience = 10L, min_delta = 0, device = "cpu",
  seed = NULL, pretrained = NULL, freeze_encoder = FALSE, concat = FALSE,
  prefix = "ae_")
```

## Arguments

- latent_dim:

  Integer, size of latent representation.

- hidden_layers:

  Integer vector of hidden layer sizes.

- activation:

  Activation function.

- dropout:

  Dropout rate between 0 and 1.

- epochs:

  Number of training epochs.

- batch_size:

  Batch size.

- lr:

  Learning rate.

- weight_decay:

  L2 weight decay.

- loss:

  Loss function: "mse" or "mae".

- validation_split:

  Fraction of data for validation.

- early_stopping:

  Logical, enable early stopping.

- patience:

  Early stopping patience.

- min_delta:

  Minimum improvement to reset patience.

- device:

  "cpu" or "cuda".

- seed:

  Optional random seed.

- pretrained:

  Optional path to a saved autoencoder state.

- freeze_encoder:

  Logical, freeze encoder weights during training.

- concat:

  Logical, if TRUE returns original + latent features.

- prefix:

  Prefix for latent feature names.

## Value

A PipeOp object.
