# Score technology-stratified rCLR

Score technology-stratified rCLR

## Usage

``` r
score_tech_stratified_rclr(
  expr,
  meta,
  panel,
  technology_col = "technology",
  fit = NULL,
  train_expr = NULL,
  train_meta = NULL,
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

- technology_col:

  Column identifying technology class.

- fit:

  Optional fit from \`fit_tech_stratified_rclr()\`.

- train_expr:

  Optional training matrix when \`fit\` is NULL.

- train_meta:

  Optional training metadata when \`fit\` is NULL.

- pseudocount:

  Optional additive pseudocount before log transform.

## Value

Numeric score vector.
