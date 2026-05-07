# Summarize Within-Provenance Case-Control Blocks

Counts cases and controls inside provenance blocks and marks whether any
block is large enough, and documented enough, to support a positive
disease-contrast claim. Missing direct-control evidence is treated
fail-closed.

## Usage

``` r
os_within_provenance_blocks(
  data,
  outcome,
  provenance_block,
  direct_control_evidence = NULL,
  positive = NULL,
  min_cases = 30L,
  min_controls = 30L
)
```

## Arguments

- data:

  Data frame containing outcome and provenance columns.

- outcome:

  Column name for the binary outcome.

- provenance_block:

  Column name defining provenance blocks.

- direct_control_evidence:

  Optional logical or character vector, or column name, indicating
  direct physical-provenance evidence for control rows.

- positive:

  Optional positive-class label. Defaults to the second factor level.

- min_cases:

  Minimum cases required inside a claim-supporting block.

- min_controls:

  Minimum controls required inside a claim-supporting block.

## Value

A data frame with per-block counts and claim-support flags.
