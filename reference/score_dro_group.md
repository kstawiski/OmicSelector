# Score the Group-DRO kit x biofluid robust discriminator

Scores each row of `X` independently with the frozen model from
[`fit_dro_group`](https://kstawiski.github.io/OmicSelector/reference/fit_dro_group.md)
by delegating to
[`score_group_dro_scorer`](https://kstawiski.github.io/OmicSelector/reference/score_group_dro_scorer.md)
on the frozen primitive fit. The primitive maps each specimen to the
frozen train-only CLR (per-row pseudocount and pivot centring), scales
it with the frozen median / MAD constants, and applies the frozen
logistic head; the returned logistic probability is already oriented so
larger = more case-like. The score of a row depends only on that row and
the frozen model.

## Usage

``` r
score_dro_group(model, X, meta = NULL)
```

## Arguments

- model:

  A `dro_group_model` object from
  [`fit_dro_group`](https://kstawiski.github.io/OmicSelector/reference/fit_dro_group.md).

- X:

  Numeric matrix (samples \\\times\\ features) of non-negative
  abundances with named feature columns.

- meta:

  Optional per-sample metadata. Forwarded to the primitive only as a
  diagnostic (it attaches test-group audit attributes) and never affects
  the numeric score, which is returned as a plain finite numeric vector.

## Value

Plain finite numeric vector of length `nrow(X)`; larger values are more
case-like.

## See also

[`fit_dro_group`](https://kstawiski.github.io/OmicSelector/reference/fit_dro_group.md),
[`score_group_dro_scorer`](https://kstawiski.github.io/OmicSelector/reference/score_group_dro_scorer.md)

## Examples

``` r
if (FALSE) { # \dontrun{
model <- fit_dro_group(X, y, meta)
score_dro_group(model, X[1, , drop = FALSE])
} # }
```
