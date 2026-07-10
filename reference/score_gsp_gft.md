# Score specimens with a frozen graph-Fourier discriminator

Scores each row of `X` independently with a frozen `gsp_gft_model`. For
each specimen, the scorer computes a per-sample rCLR vector over the
model features that are present in `X`, aligns that vector to the frozen
training feature universe with absent features set to 0, projects it
onto the frozen low-frequency graph-Fourier basis, applies the frozen
coordinate standardization, and evaluates the frozen linear head.

## Usage

``` r
score_gsp_gft(model, X, meta = NULL)
```

## Arguments

- model:

  A `gsp_gft_model` object returned by
  [`fit_gsp_gft`](https://kstawiski.github.io/OmicSelector/reference/fit_gsp_gft.md).

- X:

  Numeric matrix (samples \\\times\\ features) of non-negative abundance
  values. Columns must be named feature ids.

- meta:

  Optional per-sample metadata. Accepted for interface uniformity and
  ignored by this method.

## Value

Plain finite numeric vector of length `nrow(X)`. Larger values are more
case-like. If fewer than `model$hp$min_features` model features are
present, the documented neutral score `0` is returned for every row.

## Details

The score of row \\i\\ depends only on `X[i, ]` and the frozen model.
The scored batch contributes no centering, normalization, graph update,
bandwidth, basis, or head statistic. The rCLR representation is exactly
invariant to multiplying one specimen by a positive scalar, because
\\\log(c v_j) - mean_k \log(c v_k) = \log(v_j) - mean_k \log(v_k)\\ over
that specimen's own positive entries.

## References

Shuman DI, Narang SK, Frossard P, Ortega A, Vandergheynst P. (2013) The
emerging field of signal processing on graphs. *IEEE Signal Processing
Magazine* 30: 83-98.

Chung FRK. (1997) *Spectral Graph Theory*. American Mathematical
Society.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
X <- matrix(stats::rgamma(60 * 30, shape = 3), nrow = 60,
            dimnames = list(NULL, paste0("miR-", seq_len(30))))
y <- rep(c(0, 1), each = 30)
X[y == 1, 1:6] <- X[y == 1, 1:6] + 20
model <- fit_gsp_gft(X, y)
score_gsp_gft(model, X[1, , drop = FALSE])
} # }
```
