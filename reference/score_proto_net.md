# Score the prototypical-network single-sample discriminator

Scores each row of `X` independently with the frozen model from
[`fit_proto_net`](https://kstawiski.github.io/OmicSelector/reference/fit_proto_net.md),
in PURE R (no python). Each query is mapped to the per-sample robust CLR
over the frozen feature universe (absent universe features carry the
neutral rCLR value `0`), embedded by the exported MLP forward, and
scored by the difference of squared-Euclidean distances to the frozen
control and case prototypes. Larger = more case-like.

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
score_proto_net(model, X, meta = NULL)
```

## Arguments

- model:

  A `proto_net_model` object from
  [`fit_proto_net`](https://kstawiski.github.io/OmicSelector/reference/fit_proto_net.md).

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

Snell J, Swersky K, Zemel R. (2017) Prototypical Networks for Few-shot
Learning. *NeurIPS* 30. arXiv:1703.05175.

## Examples

``` r
if (FALSE) { # \dontrun{
model <- fit_proto_net(X, y)
score_proto_net(model, X[1, , drop = FALSE])
} # }
```
