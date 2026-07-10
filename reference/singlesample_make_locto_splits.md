# Build leave-one-cancer-type-out splits with same-block exclusion

Build leave-one-cancer-type-out splits with same-block exclusion

## Usage

``` r
singlesample_make_locto_splits(
  sample_meta,
  cancer_col = "cancer_type",
  cohort_col = "cohort",
  provenance_block_col = "provenance_block",
  min_train_cohorts = 1L
)
```

## Arguments

- sample_meta:

  Data frame with one row per sample.

- cancer_col:

  Column identifying cancer type or endpoint family.

- cohort_col:

  Column identifying cohorts.

- provenance_block_col:

  Optional column identifying specimen-provenance blocks. When present,
  samples in the held-out block are excluded from training.

- min_train_cohorts:

  Minimum number of training cohorts required for an evaluable split.

## Value

A list of split objects with train/test indices and status.
