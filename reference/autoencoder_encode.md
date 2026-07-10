# Encode Features with Autoencoder

Applies a fitted autoencoder to encode features into latent space.

## Usage

``` r
autoencoder_encode(model, x, concat = FALSE, prefix = "ae_")
```

## Arguments

- model:

  OmicAutoencoder object.

- x:

  Numeric matrix or data.frame (samples x features).

- concat:

  Logical, if TRUE returns original + latent features.

- prefix:

  Prefix for latent feature names.

## Value

A data.table with encoded features.
