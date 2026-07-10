# Build leave-one-cohort-out splits with same-block exclusion

Build leave-one-cohort-out splits with same-block exclusion

## Usage

``` r
singlesample_make_loco_splits(
  sample_meta,
  cohort_col = "cohort",
  provenance_block_col = "provenance_block",
  min_train_cohorts = 1L
)
```

## Arguments

- sample_meta:

  Data frame with one row per sample.

- cohort_col:

  Column identifying cohorts.

- provenance_block_col:

  Optional column identifying specimen-provenance blocks. When present,
  samples in the held-out block are excluded from training.

- min_train_cohorts:

  Minimum number of training cohorts required for an evaluable split.

## Value

A list of split objects with train/test indices and status.
