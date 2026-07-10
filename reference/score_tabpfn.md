# Score the TabPFN-v2 in-context discriminator (default row-by-row)

Scores each row of `X` independently with the frozen model from
[`fit_tabpfn`](https://kstawiski.github.io/OmicSelector/reference/fit_tabpfn.md).
Each query is mapped to the per-sample robust CLR over the frozen
feature universe (absent universe features carry the neutral rCLR value
`0`), then classified by TabPFN against the FROZEN training context. By
default, queries are passed to `predict_proba` ONE ROW AT A TIME (the
\\n=1\\ forward path used uniformly); the score is the case posterior
\\P(y = 1)\\ in \\\[0, 1\]\\, larger = more case-like.

Queries with fewer than `model$hp$min_features` present universe
features, or an all-zero / empty specimen, return the neutral
probability `0.5`. The score of a row depends only on that row and the
frozen model, and is exactly invariant to per-specimen positive scaling.

## Usage

``` r
score_tabpfn(model, X, meta = NULL)
```

## Arguments

- model:

  A `tabpfn_model` object from
  [`fit_tabpfn`](https://kstawiski.github.io/OmicSelector/reference/fit_tabpfn.md).

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
over balanced chunks of query rows. This path is leakage-free and
AUC-faithful (observed \\\|dAUC\| = 7.8\mathrm{e}{-5}\\ at production
scale; specimen ranks/verdicts unchanged), but is not bit-identical to
the default \\n = 1\\ row-by-row path because of a deterministic
batch-size numerical effect. The default row-by-row single-specimen
deployment path is unchanged.

## References

Hollmann N, et al. (2025) Accurate predictions on small data with a
tabular foundation model. *Nature* 637:319-326.

## Examples

``` r
if (FALSE) { # \dontrun{
model <- fit_tabpfn(X, y)
score_tabpfn(model, X[1, , drop = FALSE])
} # }
```
