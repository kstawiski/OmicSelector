# Score boosted rank trees one specimen at a time

Score boosted rank trees one specimen at a time

## Usage

``` r
score_reo_rankboost(model, X, meta = NULL)
```

## Arguments

- model:

  A \`reo_rankboost_model\` returned by \[fit_reo_rankboost()\].

- X:

  Numeric matrix of samples by the complete frozen feature universe.

- meta:

  Optional row-aligned metadata, accepted for canonical dispatch and
  otherwise ignored.

## Value

Numeric case probabilities, one per row of \`X\`.
