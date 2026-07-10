# Score curvature-scale-space (inv-css) within-sample discriminator

Scores each row of `X` independently with the frozen model learned by
[`fit_inv_css`](https://kstawiski.github.io/OmicSelector/reference/fit_inv_css.md).
The model-universe features present in `X`
(`intersect(model$feature_universe, colnames(X))`) define the present
substrate; for one specimen its within-sample rCLR profile over exactly
those present features is resampled to the frozen length, summarised by
the frozen curvature-scale-space descriptor (same scale ladder, kernels
and difference stencils as at fitting), standardised by the frozen
`center`/`scale` and passed through the frozen linear head. Larger
values are more case-like.

Scoring uses only each row's own values plus the frozen model – no
test-batch renormalisation, no cross-row coupling, no statistic
estimated from `X`. If fewer than `model$hp$min_features` model-universe
features are present, the documented neutral score `0` is returned for
every row; a single specimen whose own number of present rCLR
coordinates falls below `model$hp$min_features` likewise scores the
neutral `0`.

## Usage

``` r
score_inv_css(model, X, meta = NULL)
```

## Arguments

- model:

  An `inv_css_model` object returned by
  [`fit_inv_css`](https://kstawiski.github.io/OmicSelector/reference/fit_inv_css.md).

- X:

  Numeric matrix (samples \\\times\\ features) of non-negative abundance
  values. Columns must be named feature ids.

- meta:

  Optional per-sample metadata. Accepted for interface uniformity and
  ignored by this method.

## Value

Plain finite numeric vector of length `nrow(X)`. Larger values are more
case-like. The score of a row depends only on that row's values and the
frozen model, and is exactly invariant to per-sample positive scaling.

## References

Mokhtarian F, Mackworth AK. (1992) A theory of multiscale,
curvature-based shape representation for planar curves. *IEEE
Transactions on Pattern Analysis and Machine Intelligence* 14(8):
789-805.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
p <- 60
features <- paste0("miR-", sprintf("%03d", seq_len(p)))
n <- 160
y <- rep(c(0, 1), each = n / 2)
X <- matrix(0, n, p, dimnames = list(NULL, features))
ctrl_mu <- seq(2, 5, length.out = p)
case_mu <- ifelse(seq_len(p) <= p / 2, 2.75, 4.25)
for (i in seq_len(n)) {
  mu <- if (y[i] == 1) case_mu else ctrl_mu
  X[i, ] <- exp(stats::rnorm(p, mean = mu, sd = 0.2))
}
model <- fit_inv_css(X, y)
score_inv_css(model, X[1, , drop = FALSE])
} # }
```
