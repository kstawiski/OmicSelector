# OmicSelector_combat

Use combat to fight batch effect.

## Usage

``` r
OmicSelector_combat(
  danex,
  metadane = metadane,
  model = c("~ Batch", "~ Batch + Class")
)
```

## Arguments

- danex:

  Matrix with miRNA normalized expression values (e.g. TPM, deltaCt)
  with miRNAs in columns and cases in rows.

- metadane:

  Metadata with \`Class\` and \`Batch\` variable.

- model:

  Model of correction (application of covariates).

## Value

Batch-corrected dataset.
