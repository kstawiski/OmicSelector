# Clustered Bootstrap CI for AUC

Bootstrap the AUC of a set of predictions respecting within-cluster
dependence (e.g., duplicated specimens with the same internal ID across
accessions; repeated measures from the same patient). Rows with the same
cluster ID are resampled as a unit, which preserves the within-cluster
correlation structure and prevents anti-conservative CIs.

## Usage

``` r
os_clustered_bootstrap_auc(
  y,
  predictions,
  cluster_id,
  n_boot = 2000L,
  stratify_by = NULL,
  ci_type = c("two-sided", "one-sided-lower"),
  alpha = 0.05,
  seed = 42L
)
```

## Arguments

- y:

  Binary outcome vector.

- predictions:

  Numeric predicted probabilities / scores.

- cluster_id:

  Vector of cluster identifiers (length n).

- n_boot:

  Number of bootstrap replicates. Default 2000.

- stratify_by:

  Optional factor for stratified resampling (e.g., fold x class). Strata
  are preserved at the cluster level.

- ci_type:

  One of `"two-sided"` or `"one-sided-lower"`. Default `"two-sided"`.

- alpha:

  Significance level (default 0.05).

- seed:

  Random seed.

## Value

A list with the point AUC and bootstrap CI bounds.
