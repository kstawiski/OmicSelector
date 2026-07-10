# Fit cross-technology harmonized rCLR anchors

Fit cross-technology harmonized rCLR anchors

## Usage

``` r
fit_cross_tech_harmonized_rclr(
  train_expr,
  train_meta,
  panel,
  technology_col = "technology",
  required_technologies = c("NGS", "Toray", "qPCR", "NanoString"),
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

  Technology classes required for selection.

- min_detection_rate:

  Minimum per-technology detection rate.

- n_anchor:

  Maximum number of anchors to select.

- pseudocount:

  Optional additive pseudocount before log transform.

## Value

Fit object consumed by \`predict_cross_tech_harmonized_rclr()\`.
