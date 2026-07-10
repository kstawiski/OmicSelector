# Score the TabICL in-context discriminator (default row-by-row)

Scores each row of `X` independently with the frozen model from
[`fit_tabicl`](https://kstawiski.github.io/OmicSelector/reference/fit_tabicl.md).
Each query is mapped to the per-sample robust CLR over the frozen
feature universe (absent universe features carry the neutral rCLR value
`0`), then classified by TabICL against the FROZEN training context. By
default, queries are passed to `predict_proba` ONE ROW AT A TIME (the
\\n=1\\ forward path used uniformly); the score is the case posterior
\\P(y = 1)\\ in \\\[0, 1\]\\, larger = more case-like.

Queries with fewer than `model$hp$min_features` present universe
features, or an all-zero / empty specimen, return the neutral
probability `0.5`. The score of a row depends only on that row and the
frozen model, and is exactly invariant to per-specimen positive scaling.

## Usage

``` r
score_tabicl(model, X, meta = NULL)
```

## Arguments

- model:

  A `tabicl_model` object from
  [`fit_tabicl`](https://kstawiski.github.io/OmicSelector/reference/fit_tabicl.md).

- X:

  Numeric matrix (samples \\\times\\ features) of non-negative
  abundances with named feature columns.

- meta:

  Optional per-sample metadata. Accepted for interface uniformity and
  ignored by this method.

## Value

Plain finite numeric vector of length `nrow(X)` in \\\[0, 1\]\\; larger
values are more case-like.

## Details

A benchmark-only `score_batch` option evaluates `predict_proba` once
over balanced chunks of query rows; this is numerically identical to the
default row-by-row single-sample path because the model is
row-order-invariant and queries are conditionally independent given the
frozen context, while the \\n = 1\\ single-specimen deployment path is
unchanged.

## References

Qu J, Holzmuller D, Varoquaux G, Le Morvan M. (2025) TabICL: A Tabular
Foundation Model for In-Context Learning on Large Data. *ICML*;
arXiv:2502.05564.

## Examples

``` r
if (FALSE) { # \dontrun{
model <- fit_tabicl(X, y)
score_tabicl(model, X[1, , drop = FALSE])
} # }
```
