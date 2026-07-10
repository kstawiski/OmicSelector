# Score ordinal Local Binary Pattern (inv-olbp) within-sample discriminator

Scores each row of `X` independently with the frozen model learned by
[`fit_inv_olbp`](https://kstawiski.github.io/OmicSelector/reference/fit_inv_olbp.md).
The model-universe features present in `X` are taken in the frozen
canonical order, and for one specimen its within-sample ranks over
exactly those present features are turned into the ordinal Local Binary
Pattern histogram (same neighbour offsets, same bit packing, same
\\2^P\\ bins as at fitting; the circular neighbour wrap is recomputed on
the present length \\M\\). The histogram is standardised by the frozen
`center`/`scale` and passed through the frozen linear head. Larger
values are more case-like.

Scoring uses only each row's own values plus the frozen model – no
test-batch renormalisation, no cross-row coupling, no statistic
estimated from `X`. If the number of present model-universe features is
below `model$hp$min_features` or not greater than the neighbour count
\\P\\ (so a circular OLBP code is undefined), the documented neutral
score `0` is returned for every row.

## Usage

``` r
score_inv_olbp(model, X, meta = NULL)
```

## Arguments

- model:

  An `inv_olbp_model` object returned by
  [`fit_inv_olbp`](https://kstawiski.github.io/OmicSelector/reference/fit_inv_olbp.md).

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

Ojala T, Pietikainen M, Maenpaa T. (2002) Multiresolution gray-scale and
rotation invariant texture classification with local binary patterns.
*IEEE Transactions on Pattern Analysis and Machine Intelligence* 24(7):
971-987.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
p <- 40
X <- matrix(stats::rgamma(120 * p, shape = 20, rate = 0.2), nrow = 120,
            dimnames = list(NULL, paste0("miR-", sprintf("%03d", seq_len(p)))))
y <- rep(c(0, 1), each = 60)
X[y == 1, 11:16] <- X[y == 1, 11:16] + 150
model <- fit_inv_olbp(X, y)
score_inv_olbp(model, X[1, , drop = FALSE])
} # }
```
