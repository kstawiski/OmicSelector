# Fit biofluid-stratified rCLR centering

Selects rCLR centering features separately within each training biofluid
class. The held-out score is the signed sum of panel log-ratios centered
by the denominator selected for the sample's biofluid class.

## Usage

``` r
fit_biofluid_stratified_rclr(
  train_expr,
  train_meta,
  panel,
  biofluid_col = "biofluid",
  pseudocount = NULL,
  exclude_features = .bfa_default_blood_cell_mirnas(),
  trim_upper = 0.1,
  trim_lower = 0.05,
  min_centering_size = 8L
)
```

## Arguments

- train_expr:

  Numeric matrix, samples x features.

- train_meta:

  Data frame aligned to \`train_expr\`.

- panel:

  Character vector of panel features.

- biofluid_col:

  Column in \`train_meta\` with canonical or raw biofluid.

- pseudocount:

  Optional additive pseudocount before log transform.

- exclude_features:

  Features excluded from rCLR denominator selection.

- trim_upper, trim_lower:

  Feature trimming fractions for rCLR denominator selection.

- min_centering_size:

  Minimum denominator size before fallback.

## Value

A \`biofluid_stratified_rclr_fit\`.
