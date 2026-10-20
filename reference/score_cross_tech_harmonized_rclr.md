# Score cross-technology harmonized rCLR

Score cross-technology harmonized rCLR

## Usage

``` r
score_cross_tech_harmonized_rclr(
  expr,
  meta,
  panel,
  technology_col = "technology",
  fit = NULL,
  train_expr = NULL,
  train_meta = NULL,
  required_technologies = NULL,
  min_detection_rate = 0.05,
  n_anchor = 20L,
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

- required_technologies:

  Technology classes required for selection. \`NULL\` (default) uses the
  known classes present in training metadata after alias mapping,
  excluding \`unknown\`. Fewer than two known classes raises
  \`singlesample_input_error\` with code
  \`INSUFFICIENT_TECHNOLOGY_CLASSES\`. An explicit vector retains the
  requested classes, including absent classes (which prevent anchor
  selection); an explicit single class remains allowed.

- min_detection_rate:

  Minimum per-technology detection rate.

- n_anchor:

  Maximum number of anchors to select.

- pseudocount:

  Optional additive pseudocount before log transform.

## Value

Numeric score vector.
