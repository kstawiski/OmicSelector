# Compute Nogueira Stability Index

Computes the Nogueira Stability Index from a list of selected feature
sets across multiple cross-validation folds.

## Usage

``` r
compute_nogueira_stability(feature_sets, all_features)
```

## Arguments

- feature_sets:

  A list where each element is a character vector of selected feature
  names from one CV fold

- all_features:

  Character vector of all candidate feature names (the universe of
  features before selection)

## Value

A list containing: - nogueira_index: The stability index (0 to 1, higher
= more stable) - selection_frequency: Named vector of per-feature
selection frequencies - avg_subset_size: Average number of features
selected per fold - subset_sizes: Vector of subset sizes per fold -
consensus_features: Features selected in ALL folds - n_folds: Number of
folds analyzed - n_features: Total number of candidate features -
warning: Any warning messages (e.g., degenerate cases)

## Examples

``` r
if (FALSE) { # \dontrun{
# Feature sets from 5-fold CV
sets <- list(
  c("gene1", "gene2", "gene3"),
  c("gene1", "gene2", "gene4"),
  c("gene1", "gene2", "gene3"),
  c("gene1", "gene3", "gene5"),
  c("gene1", "gene2", "gene3")
)
all_genes <- paste0("gene", 1:100)

stability <- compute_nogueira_stability(sets, all_genes)
print(stability$nogueira_index)
} # }
```
