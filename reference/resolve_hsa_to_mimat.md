# Resolve canonical hsa-miR-\* names to MIMAT accession IDs

Resolve canonical hsa-miR-\* names to MIMAT accession IDs

## Usage

``` r
resolve_hsa_to_mimat(names)
```

## Arguments

- names:

  Character vector of canonical hsa-miR-\* mature miRNA names, e.g.
  `c("hsa-let-7a-5p", "hsa-miR-21-5p")`.

## Value

Character vector (same length as `names`) of MIMAT accession IDs.
Entries not found in the lookup table are returned as `NA_character_`.

## Details

Matching is case-insensitive on input; the returned value is always in
canonical case (MIMAT IDs uppercase, e.g. `MIMAT0000062`; hsa-miR names
lowercase, e.g. `hsa-let-7a-5p`). Unknown inputs return `NA_character_`.

## Examples

``` r
resolve_hsa_to_mimat(c("hsa-let-7a-5p", "hsa-not-real"))
#> [1] "MIMAT0000062" NA            
# [1] "MIMAT0000062" NA
```
