# Fit cross-technology harmonized rCLR anchors

Fit cross-technology harmonized rCLR anchors

## Usage

``` r
fit_cross_tech_harmonized_rclr(
  train_expr,
  train_meta,
  panel,
  technology_col = "technology",
  required_technologies = NULL,
  min_detection_rate = 0.05,
  n_anchor = 20L,
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

Fit object consumed by \`predict_cross_tech_harmonized_rclr()\`.
