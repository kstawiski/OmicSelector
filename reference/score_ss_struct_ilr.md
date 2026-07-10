# Score the structurally-informed ILR single-sample scorer (ss-struct-ilr)

Scores each row of `X` independently with the frozen
[`fit_ss_struct_ilr`](https://kstawiski.github.io/OmicSelector/reference/fit_ss_struct_ilr.md)
model. Each specimen is restricted to the frozen training feature
universe present in `X`, ILR balances are computed from the specimen's
own values via
[`ws_balance_ilr`](https://kstawiski.github.io/OmicSelector/reference/ws_balance_ilr.md)
(the per-sample pseudocount makes balances scale-invariant),
standardised with the frozen train mean/sd, and combined with the frozen
diagonal-LDA weights. No scored-batch centering, quantiles, or
renormalisation is used, so the scorer is row-equivariant
([`singlesample_assert_row_equivariant`](https://kstawiski.github.io/OmicSelector/reference/singlesample_assert_row_equivariant.md)).
Specimens with no constructible balance receive the neutral score 0.

## Usage

``` r
score_ss_struct_ilr(model, X, meta = NULL)
```

## Arguments

- model:

  A `ss_struct_ilr_model` from
  [`fit_ss_struct_ilr`](https://kstawiski.github.io/OmicSelector/reference/fit_ss_struct_ilr.md).

- X:

  Numeric matrix (samples x features) of non-negative abundances; named
  cols.

- meta:

  Optional per-sample metadata; accepted for interface uniformity,
  ignored.

## Value

Finite numeric vector of length `nrow(X)`; higher = more case-like.

## See also

[`fit_ss_struct_ilr`](https://kstawiski.github.io/OmicSelector/reference/fit_ss_struct_ilr.md)
