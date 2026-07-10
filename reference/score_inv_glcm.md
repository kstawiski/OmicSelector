# Score GLCM Haralick texture (inv-glcm) within-sample discriminator

Scores each row of `X` independently with the frozen model from
[`fit_inv_glcm`](https://kstawiski.github.io/OmicSelector/reference/fit_inv_glcm.md).
The model-universe features present in `X` are taken in the frozen
canonical order; for one specimen its within-sample rCLR landscape is
quantised with the frozen edges, a 1-D GLCM is built, the Haralick
descriptor is computed, standardised by the frozen head, and the frozen
LDA score is returned. Larger = more case-like.

Scoring uses only each row's own values plus the frozen model; no
test-batch renormalisation, no cross-row coupling, no statistic
estimated from `X`. If fewer than `model$hp$min_features` model-universe
features are present in `X`, or a single specimen has fewer than
`model$hp$min_features` positive values, the documented neutral score
`0` is returned.

## Usage

``` r
score_inv_glcm(model, X, meta = NULL)
```

## Arguments

- model:

  An `inv_glcm_model` object returned by
  [`fit_inv_glcm`](https://kstawiski.github.io/OmicSelector/reference/fit_inv_glcm.md).

- X:

  Numeric matrix (samples \\\times\\ features) of non-negative abundance
  values. Columns must be named feature ids.

- meta:

  Optional per-sample metadata. Accepted for interface uniformity;
  ignored.

## Value

Plain finite numeric vector of length `nrow(X)`. Larger values are more
case-like. Each score depends only on that row's values and the frozen
model. It is invariant to per-sample positive scaling: the rCLR
landscape is exactly scale-invariant, so the frozen-edge quantised
levels and hence the score are bit-identical, except on the measure-zero
set where an rCLR value lies within floating-point rounding of a frozen
edge.

## References

Haralick RM, Shanmugam K, Dinstein I. (1973) Textural features for image
classification. *IEEE Transactions on Systems, Man, and Cybernetics*
SMC-3(6): 610-621.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
p <- 50
X <- matrix(stats::rgamma(120 * p, shape = 20, rate = 0.2), nrow = 120,
            dimnames = list(NULL, paste0("miR-", sprintf("%03d", seq_len(p)))))
y <- rep(c(0, 1), each = 60)
X[y == 1, 11:20] <- X[y == 1, 11:20] + 150
model <- fit_inv_glcm(X, y)
score_inv_glcm(model, X[1, , drop = FALSE])
} # }
```
