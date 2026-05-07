# Rename miRNA features in a matrix or named vector

Convenience wrapper around
[`resolve_mirna_aliases`](https://kstawiski.github.io/OmicSelector/reference/resolve_mirna_aliases.md)
that accepts a samples × features matrix or a named numeric vector and
renames the features (column names of the matrix, or
[`names()`](https://rdrr.io/r/base/names.html) of the vector) to the
target namespace. Unresolved features are either dropped
(`keep_unresolved = FALSE`) or retained with their original name
(`keep_unresolved = TRUE`, default).

## Usage

``` r
apply_mirna_aliases(
  x,
  target_namespace = c("mirna_name", "mimat"),
  table = mirna_alias_table(),
  keep_unresolved = TRUE,
  verbose = FALSE,
  ...
)
```

## Arguments

- x:

  A samples × features `matrix` with column names, or a named numeric
  vector.

- target_namespace:

  One of `"mirna_name"` (default) or `"mimat"`. Passed to
  [`resolve_mirna_aliases`](https://kstawiski.github.io/OmicSelector/reference/resolve_mirna_aliases.md).

- table:

  Alias table. Defaults to
  [`mirna_alias_table()`](https://kstawiski.github.io/OmicSelector/reference/mirna_alias_table.md).

- keep_unresolved:

  Logical. If `TRUE` (default), features whose names cannot be resolved
  are kept with their original name. If `FALSE`, unresolved features are
  dropped from the output.

- verbose:

  Logical. Passed to
  [`resolve_mirna_aliases`](https://kstawiski.github.io/OmicSelector/reference/resolve_mirna_aliases.md).

- ...:

  Additional arguments passed to
  [`resolve_mirna_aliases`](https://kstawiski.github.io/OmicSelector/reference/resolve_mirna_aliases.md).

## Value

An object of the same class and structure as `x` with feature names
mapped to the target namespace. When `keep_unresolved = FALSE` and some
features are unresolved, those columns / elements are removed.

## References

Kozomara A., Birgaoanu M., Griffiths-Jones S. (2019) miRBase: from
microRNA sequences to function. *Nucleic Acids Research* 47(D1):
D155–D162. DOI:
[doi:10.1093/nar/gky1141](https://doi.org/10.1093/nar/gky1141)

## Examples

``` r
# Named numeric vector
v <- c(MIMAT0001631 = 12.3, MIMAT0000062 = 4.1, unknown_probe = 0.9)
apply_mirna_aliases(v)
#> Error in apply_mirna_aliases(v): could not find function "apply_mirna_aliases"

# Matrix
M <- matrix(runif(6), nrow = 2,
            dimnames = list(c("S1", "S2"),
                            c("MIMAT0001631", "MIMAT0000062", "junk")))
apply_mirna_aliases(M, keep_unresolved = FALSE)
#> Error in apply_mirna_aliases(M, keep_unresolved = FALSE): could not find function "apply_mirna_aliases"
```
