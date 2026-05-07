# Build Grouped Stratified Folds

Creates validation folds that keep all rows from the same biological
specimen or other grouping unit in the same fold. This is intended for
public-omics audits where duplicated specimens, repeated aliquots, or
reused biological IDs can otherwise leak across resampling splits.

## Usage

``` r
os_make_grouped_stratified_folds(
  y,
  group_id = NULL,
  strata = NULL,
  n_folds = 5L,
  seed = 1L,
  mixed_group_action = c("stop", "warn", "allow")
)
```

## Arguments

- y:

  Binary outcome vector. The second factor level is treated as positive.

- group_id:

  Optional grouping vector. If omitted, each row is its own group.

- strata:

  Optional stratification vector. Defaults to `y`.

- n_folds:

  Number of folds.

- seed:

  Random seed used for within-stratum group ordering.

- mixed_group_action:

  What to do when a group contains both outcome classes: `"stop"`,
  `"warn"`, or `"allow"`.

## Value

A list of integer test-row indices with class `"os_grouped_folds"`.
