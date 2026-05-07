# Save Autoencoder State

Saves encoder and decoder weights plus configuration to a torch file.

## Usage

``` r
autoencoder_save(model, path)
```

## Arguments

- model:

  OmicAutoencoder object.

- path:

  File path to save state (.pt recommended).

## Value

The file path (invisibly).
