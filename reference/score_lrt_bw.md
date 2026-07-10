# Score a Bures-Wasserstein Gaussian-OT class-conditional LRT discriminator

Scores each row of `X` independently with the frozen model learned by
[`fit_lrt_bw`](https://kstawiski.github.io/OmicSelector/reference/fit_lrt_bw.md).
For one specimen the abundances are mapped to the self-contained
per-sample rCLR representation \\z\\, and the score is the Gaussian
log-likelihood ratio of case vs control under the frozen BW-regularized
class Gaussians \$\$S(z) = -\tfrac{1}{2}(z-\mu\_{\mathrm{case}})^\top
\Sigma\_{\mathrm{case}}^{\mathrm{reg}\\-1}(z-\mu\_{\mathrm{case}}) -
\tfrac{1}{2}\log\det\Sigma\_{\mathrm{case}}^{\mathrm{reg}} +
\tfrac{1}{2}(z-\mu\_{\mathrm{control}})^\top
\Sigma\_{\mathrm{control}}^{\mathrm{reg}\\-1}
(z-\mu\_{\mathrm{control}}) +
\tfrac{1}{2}\log\det\Sigma\_{\mathrm{control}}^{\mathrm{reg}}.\$\$

At scoring time the present feature set is
`intersect(model$feature_universe, colnames(X))`. The class
representation (per-class mean, BW-regularized precision and log
determinant) is re-derived over exactly that present set from the frozen
raw anchors – using only frozen training data, never a scored-batch
statistic – which keeps partial-overlap scores consistent and equal to
the fit-time representation at full overlap. The score of a row
therefore depends only on that row and the frozen model, and is exactly
(to numerical tolerance) invariant to per-sample scaling. If fewer than
`model$hp$min_features` features are present, the documented neutral
score `0` is returned for every row.

## Usage

``` r
score_lrt_bw(model, X, meta = NULL)
```

## Arguments

- model:

  An `lrt_bw_model` object returned by
  [`fit_lrt_bw`](https://kstawiski.github.io/OmicSelector/reference/fit_lrt_bw.md).

- X:

  Numeric matrix (samples \\\times\\ features) of non-negative abundance
  values. Columns must be named feature ids.

- meta:

  Optional per-sample metadata. Accepted for interface uniformity and
  ignored by this method.

## Value

Plain finite numeric vector of length `nrow(X)`. Larger values are more
case-like. Scoring uses only each row's own values plus the frozen
anchors and hyperparameters.

## References

Bhatia R., Jain T., Lim Y. (2019) On the Bures-Wasserstein distance
between positive definite matrices. *Expositiones Mathematicae* 37(2):
165-191.

Takatsu A. (2011) Wasserstein geometry of Gaussian measures. *Osaka
Journal of Mathematics* 48(4): 1005-1026.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
n <- 160; p <- 30
L <- matrix(stats::rnorm(n * p, 4, 0.6), nrow = n,
            dimnames = list(NULL, paste0("miR-", seq_len(p))))
y <- rep(c(0, 1), each = n / 2)
shift <- rep(c(1.2, -1.2), length.out = 12)
L[y == 1, 1:12] <- L[y == 1, 1:12] + matrix(shift, sum(y == 1), 12, byrow = TRUE)
X <- exp(L)
model <- fit_lrt_bw(X, y)
score_lrt_bw(model, X[1, , drop = FALSE])
} # }
```
