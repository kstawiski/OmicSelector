# Score class-conditional Student-t-copula LRT discriminator

Scores each row of `X` independently with the frozen model learned by
[`fit_lrt_tcopula`](https://kstawiski.github.io/OmicSelector/reference/fit_lrt_tcopula.md).
For one specimen the abundances are mapped to the self-contained
per-sample rCLR representation \\z\\, the present features are
PIT-transformed under each class's frozen marginal and mapped to
Student-t scores \\q_c = t\_\nu^{-1}(\mathrm{PIT}\_c(z))\\, and the
score is the Student-t-copula log-density ratio \\S(z) =
\ell\_{\mathrm{case}}(q\_{\mathrm{ case}}) -
\ell\_{\mathrm{control}}(q\_{\mathrm{control}})\\ (see
[`fit_lrt_tcopula`](https://kstawiski.github.io/OmicSelector/reference/fit_lrt_tcopula.md)
for \\\ell_c\\).

At scoring time the present feature set is
`intersect(model$feature_universe, colnames(X))`. The class copula
(per-feature marginals and correlation) is re-derived over exactly that
present set from the frozen raw anchors – using only frozen training
data, never a scored-batch statistic – which keeps partial-overlap
scores consistent and equal to the fit-time representation at full
overlap. The score of a row therefore depends only on that row and the
frozen model, and is exactly invariant to per-sample scaling. If fewer
than `model$hp$min_features` features are present, the documented
neutral score `0` is returned for every row.

## Usage

``` r
score_lrt_tcopula(model, X, meta = NULL)
```

## Arguments

- model:

  An `lrt_tcopula_model` object returned by
  [`fit_lrt_tcopula`](https://kstawiski.github.io/OmicSelector/reference/fit_lrt_tcopula.md).

- X:

  Numeric matrix (samples \\\times\\ features) of non-negative abundance
  values. Columns must be named feature ids.

- meta:

  Optional per-sample metadata. Accepted for interface uniformity and
  ignored by this method.

## Value

Plain finite numeric vector of length `nrow(X)`. Larger values are more
case-like (the specimen's feature dependence matches the case copula
better than the control copula). Scoring uses only each row's own values
plus the frozen anchors and hyperparameters.

## References

Demarta S., McNeil A.J. (2005) The t copula and related copulas.
*International Statistical Review* 73(1): 111-129.

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
L[y == 1, 1:8] <- L[y == 1, 1:8] + outer(f, rep(c(1.4, -1.4), 4))
X <- exp(L)
model <- fit_lrt_tcopula(X, y, hp = list(df = 5))
score_lrt_tcopula(model, X[1, , drop = FALSE])
} # }
```
