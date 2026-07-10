# Score per-feature KDE naive-Bayes class-conditional LRT discriminator

Scores each row of `X` independently with the frozen model learned by
[`fit_lrt_nbkde`](https://kstawiski.github.io/OmicSelector/reference/fit_lrt_nbkde.md).
For one specimen the abundances are mapped to the self-contained
per-sample rCLR representation \\z\\, and the score is the naive-Bayes
class-conditional marginal log-likelihood ratio \$\$S(z) = \sum\_{j \in
\mathrm{present}} \big\[\log f\_{\mathrm{case},j}(z_j) - \log
f\_{\mathrm{control},j}(z_j)\big\],\$\$ where each \\\log f\_{c,j}\\ is
the frozen Gaussian-kernel KDE log-density of class \\c\\ for feature
\\j\\ evaluated at the specimen's own coordinate \\z_j\\ via a
numerically stable log-sum-exp.

At scoring time the present feature set is
`intersect(model$feature_universe, colnames(X))`. The per-feature KDEs
(per-class rCLR anchor values and the pooled bandwidths) are re-derived
over exactly that present set from the frozen raw anchors – using only
frozen training data, never a scored-batch statistic – which keeps
partial-overlap scores consistent and equal to the fit-time
representation at full overlap. The score of a row therefore depends
only on that row and the frozen model, and is exactly invariant to
per-sample scaling. If fewer than `model$hp$min_features` features are
present, the documented neutral score `0` is returned for every row.

## Usage

``` r
score_lrt_nbkde(model, X, meta = NULL)
```

## Arguments

- model:

  An `lrt_nbkde_model` object returned by
  [`fit_lrt_nbkde`](https://kstawiski.github.io/OmicSelector/reference/fit_lrt_nbkde.md).

- X:

  Numeric matrix (samples \\\times\\ features) of non-negative abundance
  values. Columns must be named feature ids.

- meta:

  Optional per-sample metadata. Accepted for interface uniformity and
  ignored by this method.

## Value

Plain finite numeric vector of length `nrow(X)`. Larger values are more
case-like (the specimen's per-feature marginals match the case KDEs
better than the control KDEs). Scoring uses only each row's own values
plus the frozen anchors, bandwidths and hyperparameters.

## References

Silverman B.W. (1986) *Density Estimation for Statistics and Data
Analysis*. Chapman and Hall, London.

Hand D.J., Yu K. (2001) Idiot's Bayes – not so stupid after all?
*International Statistical Review* 69(3): 385-398.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
n <- 120; p <- 30
L <- matrix(stats::rnorm(n * p, 4, 0.6), nrow = n,
            dimnames = list(NULL, paste0("miR-", seq_len(p))))
y <- rep(c(0, 1), each = n / 2)
L[y == 1, 1:8] <- L[y == 1, 1:8] + 1.2
X <- exp(L)
model <- fit_lrt_nbkde(X, y)
score_lrt_nbkde(model, X[1, , drop = FALSE])
} # }
```
