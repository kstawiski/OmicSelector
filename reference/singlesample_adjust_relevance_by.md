# Adjust the frozen directional relevance family with BY

Jointly adjusts the two prespecified directional shifted-null tests for
each frozen pair. The denominator is exactly twice the number of frozen
pairs; directions cannot be selected after observing their p-values.

## Usage

``` r
singlesample_adjust_relevance_by(
  pair_id,
  p_above_margin,
  p_below_minus_margin,
  frozen_pair_ids
)
```

## Arguments

- pair_id:

  Pair identifiers associated with the p-values.

- p_above_margin:

  Raw p-values for \`H0: effect \<= +margin\`.

- p_below_minus_margin:

  Raw p-values for \`H0: effect \>= -margin\`.

- frozen_pair_ids:

  Complete, ordered, prespecified pooled-testable pair family.

## Value

A data frame in frozen pair order containing both raw and BY-adjusted
directional p-values and the frozen directional-family size.
