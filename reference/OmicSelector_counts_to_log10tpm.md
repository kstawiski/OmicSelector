# OmicSelector_counts_to_log10tpm

Counts to log-transformed TPM-normalized counts. The funcction support
additional filter, i.e. it can leave miRNAs that appear in minimum of
\`filtr_minimalcounts\` counts in \`filtr_howmany\` samples. Usage of
the filter like that can assure that the miRNAs selected will be
detectable in qPCR.

## Usage

``` r
OmicSelector_counts_to_log10tpm(
  danex,
  metadane = metadane,
  ids = metadane$ID,
  filtr = T,
  filtr_minimalcounts = 10,
  filtr_howmany = 1/2,
  increment = 0.001
)
```

## Arguments

- danex:

  Matrix with miRNA counts with miRNAs in columns and cases in rows.

- metadane:

  Metadata with \`Class\` variable.

- ids:

  Unique identifier of samples.

- filtr:

  If expression filter should be used.

- filtr_minimalcounts:

  How many counts?

- filtr_howmany:

  In how many samples? (Please provide a percentage or proportion, e.g.
  1/2 or 1/3).

- increment:

  Increment added to TPM values before log-transformation (usually:
  0.001, so -3 will be equalt to lack of expression).

## Value

Normalized counts in \`ttpm\` format. Please note, that the function
also saves \`TPM_DGEList_filtered.rds\` to working directory, which is a
DGEList object that can be used in packages like edgeR.
