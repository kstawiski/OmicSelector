# Score specimens with a frozen single-sample selector

Score specimens with a frozen single-sample selector

## Usage

``` r
score_singlesample_selector(selector, x_new, meta = NULL)
```

## Arguments

- selector:

  A fitted \`singlesample_selector\` from
  \[fit_singlesample_selector()\].

- x_new:

  Named numeric vector, matrix, or data frame of new specimens.

- meta:

  Optional metadata, one row per new specimen.

## Value

Numeric discrimination score per new specimen. For a
\`singlesample_selector_set\`, a data frame with one column per route is
returned and an explicitly ineligible route is \`NA\`. Larger values
indicate the positive outcome class according to the training-frozen
orientation.
