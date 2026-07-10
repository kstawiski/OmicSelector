# Resolve MIMAT accession IDs to canonical hsa-miR-\* names

Resolve MIMAT accession IDs to canonical hsa-miR-\* names

## Usage

``` r
resolve_mimat_to_hsa(ids)
```

## Arguments

- ids:

  Character vector of MIMAT accession identifiers, e.g.
  `c("MIMAT0000062", "MIMAT0004481")`.

## Value

Character vector (same length as `ids`) of canonical hsa-miR-\* names.
Entries not found in the lookup table are returned as `NA_character_`.

## Details

Matching is case-insensitive on input; the returned value is always in
canonical case (MIMAT IDs uppercase, e.g. `MIMAT0000062`; hsa-miR names
lowercase, e.g. `hsa-let-7a-5p`). Unknown inputs return `NA_character_`.

## Examples

``` r
resolve_mimat_to_hsa(c("MIMAT0000062", "MIMAT9999999"))
#> [1] "hsa-let-7a-5p" NA             
# [1] "hsa-let-7a-5p" NA
```
