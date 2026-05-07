# TabDDPM Synthetic Data Generator

Generates synthetic tabular data using denoising diffusion probabilistic
models (DDPM). This is a state-of-the-art generative approach for
tabular data.

## Usage

``` r
tabddpm_generate(
  data,
  target,
  n_synthetic,
  n_steps = 1000L,
  hidden_dim = 256L,
  epochs = 1000L,
  batch_size = 32L
)
```

## Arguments

- data:

  Training data (data.frame)

- target:

  Target column name

- n_synthetic:

  Number of synthetic samples to generate

- n_steps:

  Number of diffusion steps (default: 1000)

- hidden_dim:

  Hidden layer dimension (default: 256)

- epochs:

  Training epochs (default: 1000)

- batch_size:

  Batch size (default: 32)

## Value

A data.frame of synthetic samples

## Details

TabDDPM requires the 'torch' package. The model learns the data
distribution through a denoising process and can generate realistic
synthetic samples.

## References

Kotelnikov et al. (2023). TabDDPM: Modelling Tabular Data with Diffusion
Models.

## Examples

``` r
if (FALSE) { # \dontrun{
synthetic <- tabddpm_generate(
  data = training_data,
  target = "outcome",
  n_synthetic = 100
)
} # }
```
