# Create TabTransformer Learner

Creates a TabTransformer model for tabular omics data. TabTransformer
uses self-attention to model feature interactions.

## Usage

``` r
make_tabtransformer_learner(
  n_heads = 4L,
  n_layers = 2L,
  d_model = 64L,
  dropout = 0.3,
  batch_size = 32L,
  epochs = 100L,
  learning_rate = 1e-04
)
```

## Arguments

- n_heads:

  Number of attention heads (default: 4)

- n_layers:

  Number of transformer layers (default: 2)

- d_model:

  Hidden dimension (default: 64)

- dropout:

  Dropout rate (default: 0.3)

- batch_size:

  Training batch size (default: 32)

- epochs:

  Number of epochs (default: 100)

- learning_rate:

  Learning rate (default: 0.0001)

## Value

An mlr3 Learner object

## Details

TabTransformer architecture: 1. Embed each feature into d_model
dimensions 2. Apply self-attention across features 3. Concatenate with
original features 4. MLP classification head

## References

Huang et al. (2020). TabTransformer: Tabular Data Modeling Using
Contextual Embeddings.

## Examples

``` r
if (FALSE) { # \dontrun{
learner <- make_tabtransformer_learner()
learner$train(task)
} # }
```
