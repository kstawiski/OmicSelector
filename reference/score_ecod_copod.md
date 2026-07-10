# Score the ECOD + COPOD novelty discriminator

Scores each row of `X` independently with the frozen model from
[`fit_ecod_copod`](https://kstawiski.github.io/OmicSelector/reference/fit_ecod_copod.md).
For one specimen the abundances are mapped to the per-sample robust CLR
over its own positive support; its ECOD and COPOD tail aggregates are
computed in pure R against the frozen per-feature control reference
ECDFs (never pyod); each aggregate is z-standardized with the frozen
control mean / sd, the two z-scores are averaged, and the frozen
orientation is applied so larger = more case-like.

The present feature set is
`intersect(model$feature_universe, colnames(X))`; features missing from
a specimen contribute no tail term. If fewer than
`model$hp$min_features` universe features are present, the neutral score
`0` is returned for every row. The score of a row depends only on that
row and the frozen model, and is exactly invariant to per-specimen
scaling.

## Usage

``` r
score_ecod_copod(model, X, meta = NULL)
```

## Arguments

- model:

  An `ecod_copod_model` object from
  [`fit_ecod_copod`](https://kstawiski.github.io/OmicSelector/reference/fit_ecod_copod.md).

- X:

  Numeric matrix (samples \\\times\\ features) of non-negative
  abundances with named feature columns.

- meta:

  Optional per-sample metadata. Accepted for interface uniformity and
  ignored by this method.

## Value

Plain finite numeric vector of length `nrow(X)`; larger values are more
case-like.

## Examples

``` r
if (FALSE) { # \dontrun{
model <- fit_ecod_copod(X, y)
score_ecod_copod(model, X[1, , drop = FALSE])
} # }
```
