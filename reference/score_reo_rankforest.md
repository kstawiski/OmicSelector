# Score a random rank forest one specimen at a time

Score a random rank forest one specimen at a time

## Usage

``` r
score_reo_rankforest(model, X, meta = NULL)
```

## Arguments

- model:

  A \`reo_rankforest_model\` returned by \[fit_reo_rankforest()\].

- X:

  Numeric matrix of samples by the complete frozen feature universe.

- meta:

  Optional row-aligned metadata, accepted for canonical dispatch and
  otherwise ignored.

## Value

Numeric case probabilities, one per row of \`X\`.
