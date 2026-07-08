# Compute Pareto Frontier Indices

Simple O(n^2) dominance check for Pareto optimality. Dependency-free
implementation suitable for typical candidate set sizes.

## Usage

``` r
.pareto_frontier_idx(dt, cols, maximize)
```

## Arguments

- dt:

  data.table with objective columns

- cols:

  Character vector of objective column names

- maximize:

  Logical vector (same length as cols): TRUE if higher is better

## Value

Integer vector of row indices that are Pareto-optimal
