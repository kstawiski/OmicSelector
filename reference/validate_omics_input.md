# Validate Multi-Omics Input

Validates and normalizes multi-omics input data. Accepts either a single
data.frame (backward compatible) or a named list of data matrices.

## Usage

``` r
validate_omics_input(data, target, require_common_samples = TRUE)
```

## Arguments

- data:

  Either a data.frame or a named list of data.frames/matrices

- target:

  Name of the target column

- require_common_samples:

  Logical, whether to require sample alignment

## Value

A normalized list structure with validated omics data

## Examples

``` r
if (FALSE) { # \dontrun{
# Single data.frame (backward compatible)
result <- validate_omics_input(my_data, target = "outcome")

# Multi-omics list
multi_data <- list(
  rna = rna_expression,
  mirna = mirna_expression,
  prot = protein_levels
)
result <- validate_omics_input(multi_data, target = "outcome")
} # }
```
