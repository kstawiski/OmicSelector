# Create Correlation-Based Adjacency for GNN

Builds a feature adjacency matrix from correlation structure.

## Usage

``` r
build_correlation_adjacency(
  data,
  method = "spearman",
  threshold = 0.7,
  sparse = TRUE
)
```

## Arguments

- data:

  Data matrix (samples x features)

- method:

  Correlation method ("spearman", "pearson")

- threshold:

  Minimum absolute correlation for edge (default: 0.7)

- sparse:

  Return sparse matrix (default: TRUE)

## Value

Adjacency matrix
