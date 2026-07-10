# Score biofluid-stratified rCLR

Score biofluid-stratified rCLR

## Usage

``` r
score_biofluid_stratified_rclr(
  expr,
  meta,
  panel,
  biofluid_col = "biofluid",
  fit = NULL,
  feature_weights = NULL,
  pseudocount = NULL,
  ...
)
```

## Arguments

- expr:

  Numeric matrix, samples x features.

- meta:

  Data frame aligned to \`expr\`.

- panel:

  Character vector of panel features.

- biofluid_col:

  Column in \`meta\` with canonical or raw biofluid.

- fit:

  Optional fit. If NULL, the fit is built on \`expr\` and \`meta\`.

- feature_weights:

  Optional named numeric panel weights.

- pseudocount:

  Optional additive pseudocount before log transform.

- ...:

  Additional arguments passed to \`fit_biofluid_stratified_rclr()\`.
