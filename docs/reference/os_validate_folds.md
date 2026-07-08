# Validate Grouped Folds

Validate Grouped Folds

## Usage

``` r
os_validate_folds(folds, y = NULL, group_id = NULL, strata = NULL)
```

## Arguments

- folds:

  List of integer test-row indices.

- y:

  Optional outcome vector. When supplied, coverage is checked against
  `length(y)`.

- group_id:

  Optional grouping vector used to verify that groups are not split.

- strata:

  Optional strata vector used to report fold balance.

## Value

A data frame of validation checks.
