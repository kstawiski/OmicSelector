# Score Fisher-Rao geodesic class-conditional LRT discriminator

Scores each row of `X` independently with the frozen model learned by
[`fit_ig_fisherrao`](https://kstawiski.github.io/OmicSelector/reference/fit_ig_fisherrao.md).
For one specimen \\x\\ the score is the class-conditional
Riemannian-Gaussian log-likelihood ratio \$\$S(x) = \frac{d(x,
\mu\_{ctrl})^2}{2 \sigma\_{ctrl}^2} - \frac{d(x, \mu\_{case})^2}{2
\sigma\_{case}^2},\$\$ where \\d\\ is the Fisher-Rao geodesic distance
and \\\mu_c, \sigma_c^2\\ are the class mean and dispersion. With
`model$hp$shared_dispersion = TRUE` a single pooled dispersion replaces
both \\\sigma_c^2\\, reducing the score to \\(d\_{ctrl}^2 - d\_{case}^2)
/ (2 \sigma\_{pooled}^2)\\, a Fisher-Rao nearest-centroid contrast.

At scoring time the present feature set is
`intersect(model$feature_universe, colnames(X))`. Both the specimen and
the frozen raw anchors are closed to proportions over exactly that
present set and square-root embedded (applying the same pseudocount
smoothing as at fitting), and the class means and dispersions are
re-derived over the present set from the frozen anchors – using only
frozen training data, never a scored-batch statistic – which keeps
partial-overlap distances consistent and the score exactly invariant to
per-sample scaling. The score of a row therefore depends only on that
row and the frozen model. If fewer than `model$hp$min_features` features
are present, the documented neutral score `0` is returned for every row.

## Usage

``` r
score_ig_fisherrao(model, X, meta = NULL)
```

## Arguments

- model:

  An `ig_fisherrao_model` object returned by
  [`fit_ig_fisherrao`](https://kstawiski.github.io/OmicSelector/reference/fit_ig_fisherrao.md).

- X:

  Numeric matrix (samples \\\times\\ features) of non-negative abundance
  values. Columns must be named feature ids.

- meta:

  Optional per-sample metadata. Accepted for interface uniformity and
  ignored by this method.

## Value

Plain finite numeric vector of length `nrow(X)`. Larger values are more
case-like (geometrically closer, dispersion-weighted, to the case
class). Scoring uses only each row's own values plus the frozen anchors
and hyperparameters.

## References

Rao CR. (1945) Information and accuracy attainable in the estimation of
statistical parameters. *Bulletin of the Calcutta Mathematical Society*
37: 81-91.

Amari S, Nagaoka H. (2000) *Methods of Information Geometry*. American
Mathematical Society / Oxford University Press.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
X <- matrix(stats::rgamma(60 * 30, shape = 2), nrow = 60,
            dimnames = list(NULL, paste0("miR-", seq_len(30))))
y <- rep(c(0, 1), each = 30)
X[y == 1, 1:5] <- X[y == 1, 1:5] + 20
X[y == 0, 6:10] <- X[y == 0, 6:10] + 20
model <- fit_ig_fisherrao(X, y)
score_ig_fisherrao(model, X[1, , drop = FALSE])
} # }
```
