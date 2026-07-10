# BH-FDR correction within a matched-null modality family

Applies Benjamini-Hochberg FDR correction within a single modality
family of matched-null tests. Per plan v0.3.3 section 5: matched-null
tests use BH-FDR within each of the 4 modality families separately (not
across modalities).

## Usage

``` r
singlesample_bh_fdr_correct_matched_null(results)
```

## Arguments

- results:

  List of
  [`singlesample_matched_null_benchmark`](https://kstawiski.github.io/OmicSelector/reference/singlesample_matched_null_benchmark.md)
  return objects from a single modality family.

## Value

Numeric vector of BH-FDR q-values, same length as `results`.
