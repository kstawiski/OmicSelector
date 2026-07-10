# Score rCLR after excluding blood-cell-derived miRNAs

Removes the nominated platelet/erythrocyte miRNAs from the expression
matrix before computing the rCLR denominator and before scoring the
panel.

## Usage

``` r
score_platelet_exclusion_rclr(
  expr,
  meta,
  panel,
  biofluid_col = "biofluid",
  excluded_features = .bfa_default_blood_cell_mirnas(),
  feature_weights = NULL,
  pseudocount = NULL,
  trim_upper = 0.1,
  trim_lower = 0.05,
  min_centering_size = 8L,
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

  Biofluid column retained for a common scorer signature.

- excluded_features:

  Features removed before centering and scoring.

- feature_weights:

  Optional named numeric panel weights.

- pseudocount:

  Optional additive pseudocount before log transform.

- trim_upper, trim_lower:

  Feature trimming fractions for rCLR denominator selection.

- min_centering_size:

  Minimum denominator size before fallback.

- ...:

  Reserved for scorer-interface compatibility.
