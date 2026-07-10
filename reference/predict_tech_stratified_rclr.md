# Predict technology-stratified rCLR scores

Predict technology-stratified rCLR scores

## Usage

``` r
predict_tech_stratified_rclr(
  fit,
  expr,
  meta,
  technology_col = fit$technology_col
)
```

## Arguments

- fit:

  Fit from \`fit_tech_stratified_rclr()\`.

- expr:

  Numeric matrix, samples x features.

- meta:

  Data frame aligned to \`expr\`.

- technology_col:

  Column identifying technology class.

## Value

Numeric score vector. Attribute \`fallback_samples\` counts samples
scored with the training global center because their technology was
absent from the training fit.
