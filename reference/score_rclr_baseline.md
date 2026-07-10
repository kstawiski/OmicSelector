# Baseline rCLR panel sum using the standard trimmed denominator.

Baseline rCLR panel sum using the standard trimmed denominator.

## Usage

``` r
score_rclr_baseline(
  expr,
  meta,
  panel,
  trim_upper = 0.1,
  trim_lower = 0.05,
  exclude_features = c("hsa-miR-451a", "hsa-miR-16-5p", "hsa-miR-486-5p",
    "hsa-miR-144-3p", "hsa-miR-223-3p"),
  min_centering_size = 8L,
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

- trim_upper, trim_lower:

  Feature trimming fractions for rCLR denominator selection.

- exclude_features:

  Features excluded from denominator selection.

- min_centering_size:

  Minimum denominator size before fallback.

- pseudocount:

  Optional additive pseudocount before log transform.
