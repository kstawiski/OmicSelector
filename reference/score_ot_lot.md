# Score a linearized-optimal-transport (CDT tangent) LRT discriminator

Scores each row of `X` independently with the frozen model learned by
[`fit_ot_lot`](https://kstawiski.github.io/OmicSelector/reference/fit_ot_lot.md).
For one specimen the abundances are mapped to the self-contained
per-sample rCLR representation, its present coordinates
\\z\_{\mathrm{pos}}\\ are taken, and its CDT quantile function \\Q =
\mathrm{quantile}(z\_{\mathrm{pos}}, t_k, \mathrm{type}=1)\\ (the
empirical inverse-CDF) is formed on the frozen grid. The score is the
frozen linear functional of the CDT tangent embedding \$\$S(Q) = w^\top
(Q - Q\_{\mathrm{ref}}) + b,\$\$ larger meaning more case-like.

At scoring time the present feature set is
`intersect(model$feature_universe, colnames(X))`; the per-sample rCLR is
formed over exactly that set. The reference \\Q\_{\mathrm{ref}}\\, the
grid \\t_k\\ and the head \\(w,b)\\ are frozen from training and live in
the fixed-dimension CDT tangent space, so no present-subset
re-derivation is needed and no scored-batch statistic is used. The score
of a row therefore depends only on that row and the frozen model, and is
exactly (to numerical tolerance) invariant to per-sample scaling. If
fewer than `model$hp$min_features` feature-universe columns are present,
the documented neutral score `0` is returned for every row; a specimen
whose own present rCLR coordinates number fewer than
`model$hp$min_features` (e.g. a near-empty profile) also scores the
neutral `0`.

## Usage

``` r
score_ot_lot(model, X, meta = NULL)
```

## Arguments

- model:

  An `ot_lot_model` object returned by
  [`fit_ot_lot`](https://kstawiski.github.io/OmicSelector/reference/fit_ot_lot.md).

- X:

  Numeric matrix (samples \\\times\\ features) of non-negative abundance
  values. Columns must be named feature ids.

- meta:

  Optional per-sample metadata. Accepted for interface uniformity and
  ignored by this method.

## Value

Plain finite numeric vector of length `nrow(X)`. Larger values are more
case-like. Scoring uses only each row's own values plus the frozen
reference, grid and head.

## References

Park S.R., Kolouri S., Kundu S., Rohde G.K. (2018) The cumulative
distribution transform and linear pattern classification. *Applied and
Computational Harmonic Analysis* 45(3): 616-641.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
n <- 200; p <- 40
y <- rep(c(0, 1), each = n / 2)
E <- matrix(0, n, p, dimnames = list(NULL, paste0("miR-", seq_len(p))))
E[y == 0, ] <- stats::rnorm(sum(y == 0) * p, 0, 0.30)
E[y == 1, ] <- stats::rnorm(sum(y == 1) * p, 0, 0.60)
X <- exp(6 + E)
model <- fit_ot_lot(X, y)
score_ot_lot(model, X[1, , drop = FALSE])
} # }
```
