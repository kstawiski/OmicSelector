# OmicSelector_differential_expression_ttest

The variable performes standard differential expression analysis using
unpaired t-test with BH and Bonferonni correction. It requires
\`ttpm_features\` object, which is e.g. a matrix of log-transformed
TPM-normalized miRNAs counts with miRNAs placed in \`ttpm\` design (i.e.
columns and cases placed as rows). Classess should be passed as
\`classes\` and this should be a vector of length equal to number of
rows in \`ttpm_polfiltrze\` and contain only "Case" or "Control"
labels!! The function returns the miRNAs sorted by BH-corrected p-value.

## Usage

``` r
OmicSelector_differential_expression_ttest(
  ttpm_features,
  classes,
  mode = "auto"
)
```

## Arguments

- ttpm_features:

  matrix of log-transformed TPM-normalized miRNAs counts or other
  feature matrix with miRNAs/features placed in \`ttpm\` design (i.e.
  columns and cases placed as rows)

- classes:

  vector describing label for each case. It should contain only "Case"
  and "Control" labeles!!!!

- mode:

  use 'logtpm' for log(TPM) data or 'deltact' for qPCR deltaCt values.
  This parameters sets how the fold-change is calculated. Setting it to
  "auto" will try to read settings from var_type.txt (used in docker).

## Value

Data frame with results.
