# MIMAT accession \\\leftrightarrow\\ hsa-miR name lookup

Functions to translate between MIMAT accession identifiers (e.g.
`MIMAT0000062`) and canonical *hsa-miR-\** mature miRNA names (e.g.
`hsa-let-7a-5p`) using the bundled miRBase 22.1 lookup table.

The lookup table is shipped at `inst/extdata/mimat_hsa_lookup_v22_1.tsv`
and covers all 2 656 human mature miRNAs in miRBase release 22.1. Each
MIMAT ID maps to exactly one primary name (1-to-1; no ambiguity in the
source data).
