# Fit Autoencoder

Fits a torch autoencoder and returns an OmicAutoencoder object.

## Usage

``` r
autoencoder_fit(x, latent_dim = 32L, hidden_layers = c(128L, 64L),
  activation = c("relu", "tanh", "sigmoid", "gelu"), dropout = 0.1,
  epochs = 100L, batch_size = 32L, lr = 0.001, weight_decay = 0,
  loss = c("mse", "mae"), validation_split = 0.1,
  early_stopping = TRUE, patience = 10L, min_delta = 0,
  device = "cpu", seed = NULL, pretrained = NULL,
  freeze_encoder = FALSE)
```

## Arguments

- x:

  Numeric matrix or data.frame (samples x features).

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

  Optional path to saved autoencoder state.

- freeze_encoder:

  Logical, freeze encoder weights during training.

## Value

An `OmicAutoencoder` object.
