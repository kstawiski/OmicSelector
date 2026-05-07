# Resolve miRNA identifiers to a canonical namespace

Maps a character vector of miRNA feature identifiers — in any of the
supported input namespaces — to either the canonical mature miRNA name
(`"mirna_name"`) or the primary MIMAT accession (`"mimat"`).
Case-insensitive matching is applied to MIMAT accessions (the `MIMAT`
prefix is treated case-insensitively); mature names and alias strings
are matched with trimmed whitespace but otherwise case-sensitive.

## Usage

``` r
resolve_mirna_aliases(
  features,
  target_namespace = c("mirna_name", "mimat"),
  table = mirna_alias_table(),
  keep_unresolved = TRUE,
  verbose = FALSE
)
```

## Arguments

- features:

  Character vector of feature identifiers (any mixture of MIMAT
  accessions, mature miRNA names, precursor MI accessions, or alias
  strings).

- target_namespace:

  One of `"mirna_name"` (default) or `"mimat"`; the namespace to resolve
  to.

- table:

  A `data.frame` in the format returned by
  [`mirna_alias_table`](https://kstawiski.github.io/OmicSelector/reference/mirna_alias_table.md).
  Override this to use a custom or extended alias table.

- keep_unresolved:

  Logical. If `TRUE` (default), features with no match in `table` are
  returned as `NA`. If `FALSE`, an error is raised when any feature
  cannot be resolved.

- verbose:

  Logical. If `TRUE`, emits a message listing unresolved features.
  Default `FALSE`.

## Value

Character vector of the same length as `features`, with each entry
resolved to `target_namespace` or `NA` (if unresolved and
`keep_unresolved = TRUE`).

## References

Kozomara A., Birgaoanu M., Griffiths-Jones S. (2019) miRBase: from
microRNA sequences to function. *Nucleic Acids Research* 47(D1):
D155–D162. DOI:
[doi:10.1093/nar/gky1141](https://doi.org/10.1093/nar/gky1141)

## Examples

``` r
# Toray MIMAT IDs → canonical names
resolve_mirna_aliases(c("MIMAT0001631", "MIMAT0000418"))
#> Error in resolve_mirna_aliases(c("MIMAT0001631", "MIMAT0000418")): could not find function "resolve_mirna_aliases"

# Legacy name → canonical
resolve_mirna_aliases("hsa-miR-451")
#> Error in resolve_mirna_aliases("hsa-miR-451"): could not find function "resolve_mirna_aliases"

# Reverse: canonical name → MIMAT
resolve_mirna_aliases("hsa-miR-451a", target_namespace = "mimat")
#> Error in resolve_mirna_aliases("hsa-miR-451a", target_namespace = "mimat"): could not find function "resolve_mirna_aliases"
```
