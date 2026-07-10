# Load Autoencoder State

Loads a saved autoencoder and returns a list with the model and
configuration.

## Usage

``` r
autoencoder_load(path, device = "cpu")
```

## Arguments

- path:

  Path to a saved autoencoder state.

- device:

  "cpu" or "cuda".

## Value

An `OmicAutoencoder` object.
