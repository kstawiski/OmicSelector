# Fit DeepInsight Feature Layout

Uses t-SNE on feature profiles (each feature viewed across samples) to
project features to 2D positions. Discretizes to a grid. Resolves grid
collisions by nearest-free-cell search.

Must be FIT on training data only to avoid leakage.

Requires `Rtsne` package.

## Usage

``` r
fit_deepinsight(X_train, grid_size = 14L, perplexity = NULL, seed = 42L)
```

## Arguments

- X_train:

  Numeric matrix (n_train × p features)

- grid_size:

  Integer, size of the square output grid (default: 14)

- perplexity:

  t-SNE perplexity; if NULL, auto-chosen as min(5, p-1)

- seed:

  Random seed for reproducibility

## Value

A list with:

- `positions`: p × 2 integer matrix of (row, col) per feature

- `grid_size`: the grid dimension

## References

Sharma A et al. (2019). Scientific Reports 9:11399.

## Examples

``` r
if (FALSE) { # \dontrun{
X_train <- matrix(rnorm(70), nrow = 10, ncol = 7)
layout <- fit_deepinsight(X_train, grid_size = 7, seed = 1)
} # }
```
