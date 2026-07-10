# Score the uLSIF density-ratio LRT discriminator

Scores each row of `X` independently with the frozen model learned by
[`fit_dre_ulsif`](https://kstawiski.github.io/OmicSelector/reference/fit_dre_ulsif.md).
For one specimen the abundances are mapped to the self-contained
per-sample rCLR representation \\z\\, and the score is the log estimated
density ratio \$\$S(z) = \log\max\\\Bigl(\sum_l \alpha_l
\exp\bigl(-\lVert z - c_l\rVert^2 / (2\sigma^2)\bigr),\\
\varepsilon\Bigr),\$\$ with larger values more case-like.

At scoring time the present feature set is
`intersect(model$feature_universe, colnames(X))`. The kernel centres,
bandwidth and coefficients are re-derived over exactly that present set
from the frozen raw anchors / centres – using only frozen training data,
never a scored-batch statistic – which keeps partial-overlap scores
consistent and equal to the fit-time representation at full overlap. The
score of a row therefore depends only on that row and the frozen model,
and is exactly invariant to per-sample scaling. If fewer than
`model$hp$min_features` features are present, the documented neutral
score `0` is returned for every row.

## Usage

``` r
score_dre_ulsif(model, X, meta = NULL)
```

## Arguments

- model:

  A `dre_ulsif_model` object returned by
  [`fit_dre_ulsif`](https://kstawiski.github.io/OmicSelector/reference/fit_dre_ulsif.md).

- X:

  Numeric matrix (samples \\\times\\ features) of non-negative abundance
  values. Columns must be named feature ids.

- meta:

  Optional per-sample metadata. Accepted for interface uniformity and
  ignored by this method.

## Value

Plain finite numeric vector of length `nrow(X)`. Larger values are more
case-like (the specimen lies where the case density dominates the
control density). Scoring uses only each row's own values plus the
frozen centres, bandwidth and coefficients.

## References

Kanamori T, Hido S, Sugiyama M. (2009) A least-squares approach to
direct importance estimation. *Journal of Machine Learning Research* 10:
1391-1445.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
n <- 120; p <- 30; k <- 10
L <- matrix(stats::rnorm(n * p, 4, 0.5), nrow = n,
            dimnames = list(NULL, paste0("miR-", seq_len(p))))
y <- rep(c(0, 1), each = n / 2)
L[y == 1, seq_len(k)] <- L[y == 1, seq_len(k)] + 1.2
X <- exp(L)
model <- fit_dre_ulsif(X, y)
score_dre_ulsif(model, X[1, , drop = FALSE])
} # }
```
