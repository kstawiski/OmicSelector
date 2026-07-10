# Score class-conditional vine-copula LRT discriminator

Scores each row of `X` independently with the frozen model learned by
[`fit_lrt_vinecopula`](https://kstawiski.github.io/OmicSelector/reference/fit_lrt_vinecopula.md).
For one specimen the abundances are mapped to the self-contained
per-sample rCLR representation \\z\\, the present features are
PIT-transformed under each class's frozen marginal to
pseudo-observations \\u_c = \mathrm{PIT}\_c(z)\\, and the score is the
vine-copula log-density ratio \$\$S(z) = \log
c\_{\mathrm{case}}(u\_{\mathrm{case}}) - \log
c\_{\mathrm{control}}(u\_{\mathrm{control}}).\$\$

At scoring time the present feature set is
`intersect(model$feature_universe, colnames(X))`. The class vine
(per-feature marginals and pair-copulas) is re-derived over exactly that
present set from the frozen raw anchors – using only frozen training
data, never a scored-batch statistic – which keeps partial-overlap
scores consistent and equal to the fit-time representation at full
overlap. The score of a row therefore depends only on that row and the
frozen model, and is exactly invariant to per-sample scaling. If fewer
than `model$hp$min_features` features are present, the documented
neutral score `0` is returned for every row.

## Usage

``` r
score_lrt_vinecopula(model, X, meta = NULL)
```

## Arguments

- model:

  An `lrt_vinecopula_model` object returned by
  [`fit_lrt_vinecopula`](https://kstawiski.github.io/OmicSelector/reference/fit_lrt_vinecopula.md).

- X:

  Numeric matrix (samples \\\times\\ features) of non-negative abundance
  values. Columns must be named feature ids.

- meta:

  Optional per-sample metadata. Accepted for interface uniformity and
  ignored by this method.

## Value

Plain finite numeric vector of length `nrow(X)`. Larger values are more
case-like (the specimen's feature dependence matches the case vine
better than the control vine). Scoring uses only each row's own values
plus the frozen anchors and hyperparameters.

## References

Aas K., Czado C., Frigessi A., Bakken H. (2009) Pair-copula
constructions of multiple dependence. *Insurance: Mathematics and
Economics* 44(2): 182-198.

Joe H. (2014) *Dependence Modeling with Copulas*. Chapman and Hall/CRC.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
n <- 120; p <- 30
L <- matrix(stats::rnorm(n * p, 4, 0.6), nrow = n,
            dimnames = list(NULL, paste0("miR-", seq_len(p))))
y <- rep(c(0, 1), each = n / 2)
f <- stats::rnorm(sum(y == 1))
L[y == 1, 1:8] <- L[y == 1, 1:8] + outer(f, rep(1.4, 8))
X <- exp(L)
model <- fit_lrt_vinecopula(X, y)
score_lrt_vinecopula(model, X[1, , drop = FALSE])
} # }
```
