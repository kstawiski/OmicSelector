# Score topological persistence-image (tda-ph) within-sample discriminator

Scores each row of `X` independently with the frozen model learned by
[`fit_tda_ph`](https://kstawiski.github.io/OmicSelector/reference/fit_tda_ph.md).
The model-universe features present in `X` are taken in the frozen
canonical order, a per-specimen rCLR landscape is built over that row's
own positive support, finite 0-dimensional lower-star classes are
rasterised through the frozen exact-integral persistence-image grid, and
the resulting vector is passed through the frozen standardisation and
ridge-LDA head. Larger values are more case-like.

Scoring uses only each row's own values plus the frozen model – no
test-batch renormalisation, no cross-row coupling, no statistic
estimated from `X`. If fewer than `model$hp$min_features` model-universe
features are present, the documented neutral score `0` is returned for
every row; a specimen whose own positive-support count is below
`min_features` likewise scores the neutral `0`. A specimen with
*sufficient* support but an empty finite diagram (e.g. a strictly
monotone landscape) is *not* abstained: its all-zero persistence image
is standardised and scored by the frozen head, yielding the constant
empty-diagram baseline described in
[`fit_tda_ph`](https://kstawiski.github.io/OmicSelector/reference/fit_tda_ph.md)
(which is generally not the bare intercept \\b\\).

## Usage

``` r
score_tda_ph(model, X, meta = NULL)
```

## Arguments

- model:

  A `tda_ph_model` object returned by
  [`fit_tda_ph`](https://kstawiski.github.io/OmicSelector/reference/fit_tda_ph.md).

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

Edelsbrunner H, Letscher D, Zomorodian A. (2002) Topological persistence
and simplification. *Discrete & Computational Geometry* 28(4): 511-533.

Adams H, Emerson T, Kirby M, et al. (2017) Persistence images: a stable
vector representation of persistent homology. *Journal of Machine
Learning Research* 18(8): 1-35.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
X <- matrix(stats::rgamma(80 * 40, shape = 10, rate = 1), nrow = 80,
            dimnames = list(NULL, paste0("miR-", seq_len(40))))
y <- rep(c(0, 1), each = 40)
X[y == 1, 15:24] <- X[y == 1, 15:24] * rep(c(1.8, 0.55), 5)
model <- fit_tda_ph(X, y)
score_tda_ph(model, X[1, , drop = FALSE])
} # }
```
