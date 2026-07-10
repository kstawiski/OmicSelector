# Score the SCARF self-supervised contrastive single-sample discriminator

Scores each row of `X` independently with the frozen model from
[`fit_ai_scarf`](https://kstawiski.github.io/OmicSelector/reference/fit_ai_scarf.md),
in PURE R (no python). Each query is mapped to the per-sample robust CLR
over the frozen feature universe (absent universe features carry the
neutral rCLR value `0`), embedded by the exported encoder forward (a
CLEAN forward – NO corruption; the corruption marginal is not consulted
at score), and scored by the frozen linear-probe logit \\w \cdot z +
b\\. Larger = more case-like.

Queries with fewer than `model$hp$min_features` universe features
present (a column-overlap floor), or empty positive support over the
frozen universe on the (pre-rCLR) aligned abundances
(`!any(X_use[i, ] > 0)`), return the neutral score `0`. A FLAT
all-equal-positive composition maps to the rCLR origin but is a VALID
specimen and is scored normally. The score of a row depends only on that
row and the frozen model and is exactly invariant to per-specimen
positive scaling.

## Usage

``` r
score_ai_scarf(model, X, meta = NULL)
```

## Arguments

- model:

  An `ai_scarf_model` object from
  [`fit_ai_scarf`](https://kstawiski.github.io/OmicSelector/reference/fit_ai_scarf.md).

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

Bahri D, Jiang H, Tay Y, Metzler D. (2022) SCARF: Self-Supervised
Contrastive Learning using Random Feature Corruption. *ICLR*.
arXiv:2106.15147.

## Examples

``` r
if (FALSE) { # \dontrun{
model <- fit_ai_scarf(X, y)
score_ai_scarf(model, X[1, , drop = FALSE])
} # }
```
