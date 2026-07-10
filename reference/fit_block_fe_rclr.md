# Fit provenance-block rCLR centering moments

Computes one robust log-geometric reference vector per provenance block
from a training pool. Prediction subtracts the block-specific reference
from each sample's log abundance and sums the panel residuals.

## Usage

``` r
fit_block_fe_rclr(
  expr,
  meta,
  panel,
  block_col = "provenance_block",
  pseudocount = NULL
)
```

## Arguments

- expr:

  Numeric matrix, samples x features, training pool.

- meta:

  Data frame aligned to \`expr\`.

- panel:

  Character vector of panel features.

- block_col:

  Column in \`meta\` identifying provenance blocks.

- pseudocount:

  Optional additive pseudocount before log transform.

## Value

A fit object consumed by \`predict_block_fe_rclr()\`.
