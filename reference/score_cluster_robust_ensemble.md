# Score a cluster-robust ensemble of rCLR, ALR, and ILR-style components

The sample-level score is the row mean of robustly scaled rCLR, ALR, and
simple ILR component scores. Cluster-robust confidence intervals are not
a property of a single score vector; the benchmark script estimates them
by resampling cohort/provenance clusters across held-out folds.

## Usage

``` r
score_cluster_robust_ensemble(
  expr,
  meta,
  panel,
  cohort_col = "cohort",
  block_col = "provenance_block",
  pseudocount = NULL
)
```

## Arguments

- expr:

  Numeric matrix, samples x features.

- meta:

  Data frame aligned to \`expr\`.

- panel:

  Character vector of panel features.

- cohort_col:

  Column in \`meta\` identifying cohorts.

- block_col:

  Column in \`meta\` identifying provenance blocks.

- pseudocount:

  Optional additive pseudocount before log transform.

## Value

Numeric score vector. Attribute \`component_scores\` contains the finite
component matrix.
