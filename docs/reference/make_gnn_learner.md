# Create GNN Learner for Pathway-Aware Classification

Creates a Graph Neural Network that incorporates biological pathway
information for omics classification.

## Usage

``` r
make_gnn_learner(
  adjacency_matrix = NULL,
  n_hidden = 64L,
  n_layers = 2L,
  dropout = 0.3,
  aggregation = "mean",
  epochs = 100L
)
```

## Arguments

- adjacency_matrix:

  Adjacency matrix defining feature relationships (e.g., from pathway
  databases, correlation, or PPI networks)

- n_hidden:

  Number of hidden units (default: 64)

- n_layers:

  Number of GNN layers (default: 2)

- dropout:

  Dropout rate (default: 0.3)

- aggregation:

  Message aggregation: "mean", "sum", "max"

- epochs:

  Number of epochs (default: 100)

## Value

An mlr3 Learner object or a custom GNN wrapper

## Details

The GNN uses: - Graph Convolutional layers to propagate information -
Feature relationships from biological networks - Node-level predictions
aggregated to graph-level

## References

Kipf & Welling (2017). Semi-Supervised Classification with Graph
Convolutional Networks.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create adjacency from correlation
cor_mat <- cor(data)
adj <- (abs(cor_mat) > 0.7) * 1
diag(adj) <- 0

learner <- make_gnn_learner(adj)
learner$train(task)
} # }
```
