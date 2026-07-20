# Specimen-Duplication Detection Across Cohorts

Detects putative specimen duplicates across multiple cohorts/datasets by
either (a) exact feature equality on a chosen fingerprint matrix, or (b)
an explicit crosswalk of internal specimen IDs.

In circulating-miRNA cancer meta-analyses this is critical because the
same biological samples are often submitted to GEO under different
series accessions with different GSM IDs. Failing to dedup leads to
specimen-level train/test leakage and anti-conservative bootstrap CIs
(documented in GSE106817/GSE113486/GSE113740 for ovarian cases).

## Usage

``` r
os_detect_cross_cohort_duplicates(
  X,
  cohort,
  sample_id,
  specimen_id = NULL,
  tolerance = 0,
  min_features = NULL
)
```

## Arguments

- X:

  Numeric matrix with one row per sample, columns are features.

- cohort:

  Vector of cohort identity per sample.

- sample_id:

  Vector of sample IDs.

- specimen_id:

  Optional explicit internal specimen IDs (preferred).

- tolerance:

  Numeric tolerance for exact equality on feature vectors. Default 0
  (numeric equality without rounding).

- min_features:

  Minimum number of jointly finite features required to declare a match.
  All jointly finite features must agree within \`tolerance\`. Default:
  all columns, which therefore requires complete profiles.

## Value

A data frame with columns: `sample_id_a`, `cohort_a`, `sample_id_b`,
`cohort_b`, `specimen_id` (if provided), `match_type`.
