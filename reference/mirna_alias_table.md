# Curated miRNA alias lookup table (miRBase v22.1)

Returns a `data.frame` with one row per canonical mature miRNA and
columns for the primary mature miRNA name, MIMAT accession, optional
precursor MI accession, and a semicolon-separated list of known
alternate identifiers. The table covers the ~92 high-priority
circulating miRNAs used in the OmicSelector single-sample scoring bank:
all members of
[`ws_default_sbp`](https://kstawiski.github.io/OmicSelector/reference/ws_default_sbp.md)
and
[`ws_default_pivot_pool`](https://kstawiski.github.io/OmicSelector/reference/ws_default_pivot_pool.md),
the five canonical haemolysis markers, the Mitchell 2008 / miRBiT panel
members, and frequently-deposited Toray / FirePlex identifiers.

## Usage

``` r
mirna_alias_table()
```

## Value

A `data.frame` with columns:

- mirna_name:

  Canonical mature miRNA name (miRBase v22.1), e.g. `"hsa-miR-451a"`.

- mimat:

  Primary MIMAT accession, e.g. `"MIMAT0001631"`.

- mimat_pre:

  Precursor MI accession (character; `NA` when not unambiguously
  resolvable to a single mature arm).

- aliases:

  Semicolon-separated list of alternate identifiers including legacy
  names, abbreviated forms, and platform-specific probe name stems.

## Details

MIMAT accessions are taken from miRBase v22.1 (released March 2018).
Aliases include: legacy names from earlier miRBase releases, Affy probe
ID stems, Agilent probe name stems, and common shortened forms
encountered in published GEO depositions. This table is frozen at
package build time; downstream callers that need fresher annotations
should supply a custom table via the `table` argument of
[`resolve_mirna_aliases`](https://kstawiski.github.io/OmicSelector/reference/resolve_mirna_aliases.md).

## References

Kozomara A., Birgaoanu M., Griffiths-Jones S. (2019) miRBase: from
microRNA sequences to function. *Nucleic Acids Research* 47(D1):
D155–D162. DOI:
[doi:10.1093/nar/gky1141](https://doi.org/10.1093/nar/gky1141)

## Examples

``` r
tbl <- mirna_alias_table()
nrow(tbl)                                # number of curated miRNAs
#> [1] 88
tbl[tbl$mimat == "MIMAT0001631", ]       # hsa-miR-451a row
#>     mirna_name        mimat mimat_pre                                  aliases
#> 1 hsa-miR-451a MIMAT0001631 MI0001729 miR-451;hsa-miR-451;miR451a;MIMAT0001631
```
