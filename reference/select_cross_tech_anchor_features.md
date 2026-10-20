# Select cross-technology anchor features on a training pool

A selected anchor must be detected in every required technology class.
Selection is based only on \`train_expr\` and \`train_meta\`; held-out
samples should never be passed here.

## Usage

``` r
select_cross_tech_anchor_features(
  train_expr,
  train_meta,
  technology_col = "technology",
  required_technologies = NULL,
  min_detection_rate = 0.05,
  n_anchor = 20L
)
```

## Arguments

- train_expr:

  Numeric matrix, samples x features, training pool.

- train_meta:

  Data frame aligned to \`train_expr\`.

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

## Value

Data frame with detection rates and \`selected_yes\`.
