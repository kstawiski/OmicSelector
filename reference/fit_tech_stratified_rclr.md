# Fit technology-stratified rCLR centering moments

Fit technology-stratified rCLR centering moments

## Usage

``` r
fit_tech_stratified_rclr(
  train_expr,
  train_meta,
  panel,
  technology_col = "technology",
  pseudocount = NULL
)
```

## Arguments

- train_expr:

  Numeric matrix, samples x features, training pool.

- train_meta:

  Data frame aligned to \`train_expr\`.

- panel:

  Character vector of panel features.

- technology_col:

  Column identifying technology class.

- pseudocount:

  Optional additive pseudocount before log transform.

## Value

Fit object consumed by \`predict_tech_stratified_rclr()\`.
