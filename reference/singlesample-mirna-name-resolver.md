# Single-sample miRNA name / alias resolver (Module A, P1 — cross-platform)

Platform-portable alias resolver for circulating miRNA feature vectors.
Microarray platforms deposited in public repositories use a variety of
identifier namespaces: Toray 3D-Gene probe IDs are raw MIMAT accessions
(e.g. `MIMAT0001631`); Affymetrix miRNA-3_0 / -4_0 arrays use
Affy-internal probe IDs; Agilent miRNA arrays use Agilent probe IDs;
NanoString and FirePlex panels typically use mature-name strings of
varying vintage (miRBase v19–v22). Because
[`ws_balance_ilr`](https://kstawiski.github.io/OmicSelector/reference/ws_balance_ilr.md)
and
[`ws_alr_pivot`](https://kstawiski.github.io/OmicSelector/reference/ws_alr_pivot.md)
key the biology-frozen partition dictionary on canonical mature miRNA
names (e.g. `hsa-miR-451a`, `hsa-let-7a-5p`), a feature vector arriving
with MIMAT IDs will produce NA balances and fall to an AUC of 0.500.
This module resolves that mismatch for the ~60–100 highest-priority
circulating-miRNA features used in the single-sample scoring bank.

Three exported functions are provided:

- [`mirna_alias_table`](https://kstawiski.github.io/OmicSelector/reference/mirna_alias_table.md)
  — curated miRBase v22.1 lookup table mapping canonical mature names to
  MIMAT accessions and common aliases.

- [`resolve_mirna_aliases`](https://kstawiski.github.io/OmicSelector/reference/resolve_mirna_aliases.md)
  — maps a character vector of identifiers in any supported namespace to
  the target namespace.

- [`apply_mirna_aliases`](https://kstawiski.github.io/OmicSelector/reference/apply_mirna_aliases.md)
  — convenience wrapper that renames the columns of a samples × features
  matrix (or names of a numeric vector).

Supported input namespaces (detected automatically):

- Canonical mature miRNA name (`hsa-miR-451a`, `hsa-let-7a-5p`)

- MIMAT accession (`MIMAT0001631`) — case-insensitive prefix match

- Precursor MI accession (`MI0001729`) — resolves to the primary arm

- Legacy / alias names stored in the `aliases` column
  (semicolon-separated)

## References

Kozomara A., Birgaoanu M., Griffiths-Jones S. (2019) miRBase: from
microRNA sequences to function. *Nucleic Acids Research* 47(D1):
D155–D162. DOI:
[doi:10.1093/nar/gky1141](https://doi.org/10.1093/nar/gky1141)

Mitchell P. S., Parkin R. K., Kroh E. M., et al. (2008) Circulating
microRNAs as stable blood-based markers for cancer detection.
*Proceedings of the National Academy of Sciences USA* 105(30):
10513–10518.
