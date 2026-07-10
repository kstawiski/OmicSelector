# Predict cross-technology harmonized rCLR scores

Predict cross-technology harmonized rCLR scores

## Usage

``` r
predict_cross_tech_harmonized_rclr(
  fit,
  expr,
  meta,
  technology_col = fit$technology_col,
  min_anchors_present = 3L
)
```

## Arguments

- fit:

  Fit from \`fit_cross_tech_harmonized_rclr()\`.

- expr:

  Numeric matrix, samples x features.

- meta:

  Data frame aligned to \`expr\`; technology labels are validated but
  not used to select anchors.

- technology_col:

  Column identifying technology class.

- min_anchors_present:

  Minimum anchors required in \`expr\`.

## Value

Numeric score vector with attribute \`n_anchors_present\`.
