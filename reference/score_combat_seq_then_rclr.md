# Optional ComBat-seq followed by rCLR scoring.

Optional ComBat-seq followed by rCLR scoring.

## Usage

``` r
score_combat_seq_then_rclr(
  expr,
  meta,
  panel,
  batch_col = NULL,
  disease_col = NULL,
  allow_transductive = FALSE,
  pseudocount = 0.5
)
```

## Arguments

- expr:

  Numeric count matrix, samples x features.

- meta:

  Data frame aligned to \`expr\`.

- panel:

  Character vector of panel features.

- batch_col:

  Optional batch column. Defaults to inferred kit column.

- disease_col:

  Optional disease-label column passed to ComBat-seq.

- allow_transductive:

  Logical. Must be TRUE because this comparator is transductive and
  should not be treated as a frozen deployable scorer.

- pseudocount:

  Additive pseudocount before rCLR scoring.
