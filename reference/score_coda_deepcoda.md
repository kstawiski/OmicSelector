# Score the DeepCoDA zero-sum log-contrast single-sample discriminator

Scores each row of `X` independently with the frozen model from
[`fit_coda_deepcoda`](https://kstawiski.github.io/OmicSelector/reference/fit_coda_deepcoda.md),
in PURE R (no python). Each query is mapped to the per-sample robust CLR
over the frozen feature universe (absent universe features carry the
neutral rCLR value `0`), projected through the frozen zero-sum
log-contrast bottleneck, passed through the self-explaining
\\\theta\\-net, and scored by the DeepCoDA logit \\\sum_k
\theta_k(x)\\b_k(x) + c\\. Larger = more case-like.

Queries with fewer than `model$hp$min_features` universe features
present (a column-overlap floor), or empty positive support over the
frozen universe on the (pre-rCLR) aligned abundances
(`!any(X_use[i, ] > 0)`), return the neutral score `0`. A FLAT
all-equal-positive composition maps to the rCLR origin but is a VALID
specimen and is scored normally (its bottlenecks are all 0, so the logit
reduces to the learned bias – a genuine computed value, not a floored
0). The score of a row depends only on that row and the frozen model and
is exactly invariant to per-specimen positive scaling.

## Usage

``` r
score_coda_deepcoda(model, X, meta = NULL)
```

## Arguments

- model:

  A `coda_deepcoda_model` object from
  [`fit_coda_deepcoda`](https://kstawiski.github.io/OmicSelector/reference/fit_coda_deepcoda.md).

- X:

  Numeric matrix (samples \\\times\\ features) of non-negative
  abundances with named feature columns.

- meta:

  Optional per-sample metadata. Accepted for interface uniformity and
  ignored by this method.

## Value

Plain finite numeric vector of length `nrow(X)`; larger values are more
case-like.

## References

Quinn TP, Nguyen D, Rana S, Gupta S, Venkatesh S. (2020) DeepCoDA:
personalized interpretability for compositional health data. *ICML*,
PMLR 119. arXiv:2006.01392.

## Examples

``` r
if (FALSE) { # \dontrun{
model <- fit_coda_deepcoda(X, y)
score_coda_deepcoda(model, X[1, , drop = FALSE])
} # }
```
