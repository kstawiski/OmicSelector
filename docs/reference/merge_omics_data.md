# Merge Multi-Omics Data for Analysis

Merges multiple omics modalities into a single data.frame with
namespaced features for downstream analysis.

## Usage

``` r
merge_omics_data(omics_input, sample_subset = NULL)
```

## Arguments

- omics_input:

  A validated OmicsInput object

- sample_subset:

  Optional character vector of sample IDs to include

## Value

A merged data.frame with namespaced feature columns
