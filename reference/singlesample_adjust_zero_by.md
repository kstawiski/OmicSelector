# Adjust the frozen zero-difference test family with BY

Applies Benjamini–Yekutieli adjustment using the complete frozen pair
family as the denominator. Exactly one row per frozen pair is required;
an unavailable p-value must be represented explicitly as \`NA\` and is
retained in the output.

## Usage

``` r
singlesample_adjust_zero_by(pair_id, p_zero, frozen_pair_ids)
```

## Arguments

- pair_id:

  Pair identifiers associated with \`p_zero\`.

- p_zero:

  Raw two-sided zero-difference p-values.

- frozen_pair_ids:

  Complete, ordered, prespecified pooled-testable pair family.

## Value

A data frame in frozen family order with raw and BY-adjusted p-values
and the frozen family size.
