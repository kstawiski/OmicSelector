# Enumerate frozen within-method pairs

Enumerates every unordered pair of methods whose roster estimand is
\`"within"\`. Pair direction follows roster row order and never depends
on an outcome, score, eligibility, or performance column.

## Usage

``` r
singlesample_within_method_pairs(roster = singlesample_method_roster())
```

## Arguments

- roster:

  Method roster containing unique \`method_id\` and \`estimand\`
  columns. Defaults to \[singlesample_method_roster()\].

## Value

A data frame with one row per unordered pair and columns \`pair_id\`,
\`method_a\`, \`method_b\`, \`roster_order_a\`, and \`roster_order_b\`.
