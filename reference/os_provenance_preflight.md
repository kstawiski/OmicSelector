# Paper 3 specimen-overlap provenance pre-flight gate

Audits a candidate accession list against a manifest of known specimen-
overlap pairs and returns one of three exit conditions:

- `NO_OVERLAP`:

  No accession in the input list overlaps with any other accession in
  the input list according to the manifest, and every supplied accession
  appears (as either side of a pair) somewhere in the manifest.

- `KNOWN_OVERLAP`:

  At least one pair of input accessions is listed in the manifest. The
  cross-cohort analysis must apply block-aware multiplicity correction
  ([`singlesample_bh_fdr_correct_blocked`](https://kstawiski.github.io/OmicSelector/reference/singlesample_bh_fdr_correct_blocked.md)).

- `UNKNOWN_ACCESSION`:

  At least one input accession is not in the manifest. Audit required
  before combining with anything else. Takes precedence over
  `KNOWN_OVERLAP`.

Intended as the first stage in the Paper-3 pipeline (Methods section
"Provenance auditing and specimen-overlap handling").

## Usage

``` r
os_provenance_preflight(accessions, manifest_path = NULL)
```

## Arguments

- accessions:

  Character vector of GEO/ArrayExpress accessions (e.g.,
  `c("GSE211692", "GSE106817")`).

- manifest_path:

  Path to a TSV file with columns `acc1, acc2, block_id` (and optionally
  `fraction`, `evidence_source`). If `NULL`, uses the bundled manifest
  at
  `system.file("extdata", "provenance_manifest.tsv", package = "OmicSelector")`.

## Value

List with components:

- `status`:

  One of `"NO_OVERLAP"`, `"KNOWN_OVERLAP"`, `"UNKNOWN_ACCESSION"`.

- `overlapping_pairs`:

  `data.frame` subset of the manifest with rows whose `acc1` and `acc2`
  are both in the input accession list.

- `unknown_accessions`:

  Character vector of input accessions not present in the manifest.

- `manifest_path`:

  Resolved path of the manifest used.

## Examples

``` r
res <- os_provenance_preflight(c("GSE211692", "GSE106817"))
res$status
#> [1] "KNOWN_OVERLAP"
```
