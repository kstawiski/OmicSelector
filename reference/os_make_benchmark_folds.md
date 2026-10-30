# Reproduce the archived benchmark grouped folds

Groups are enumerated in first-row order, classified by rounded mean
outcome, shuffled within outcome (0 then 1), and assigned round-robin.

## Usage

``` r
os_make_benchmark_folds(y, group_id, k = 5L, seed = 42L)
```

## Arguments

- y:

  Numeric binary outcome vector (0 and 1).

- group_id:

  Nonmissing group identifiers, one per row.

- k:

  Number of folds, at least two.

- seed:

  Random seed.

## Value

An unnamed list of sorted integer test-row indices.

## Details

This reproduces the benchmark harness algorithm and changes the global
RNG state through `set.seed(seed)`. A mixed-outcome group with mean 0.5
is assigned to class 0. Supply original row order and the same `RNGkind`
settings (the archive uses the R default). The existing
[`os_make_grouped_stratified_folds`](https://kstawiski.github.io/OmicSelector/reference/os_make_grouped_stratified_folds.md)
balances row counts and retains its different partitions and mixed-group
checks.
