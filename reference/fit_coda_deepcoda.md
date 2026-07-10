# Fit the DeepCoDA zero-sum log-contrast single-sample discriminator

Trains a zero-sum log-contrast bottleneck (DeepCoDA eq. 1) plus a
self-explaining per-sample predictor end-to-end on the per-sample-rCLR,
universe-aligned TRAINING matrix by binary cross-entropy (via
reticulate-python torch, Adam, `BCEWithLogitsLoss`), EXPORTS the centred
(zero-sum) bottleneck weight, the \\\theta\\-net weights, and the scalar
bias to R as plain numeric matrices/vectors, and DISCARDS the python
module. The fitted model holds no external pointer; scoring is pure R.

## Usage

``` r
fit_coda_deepcoda(X_train, y_train, meta_train = NULL, hp = list())
```

## Arguments

- X_train:

  Numeric matrix (samples \\\times\\ features) of non-negative
  abundances with unique, non-empty feature names (colnames).

- y_train:

  Numeric / integer 0/1 labels (1 = case), length `nrow(X_train)`, with
  at least one case and one control.

- meta_train:

  Optional per-sample metadata. Accepted for interface uniformity and
  otherwise ignored.

- hp:

  Optional list of hyperparameters. Allowed fields: `bottleneck_dim`
  (number of zero-sum log-contrast bottlenecks \\B\\, positive integer;
  default `8L`), `theta_hidden` (positive-integer vector of
  self-explaining \\\theta\\-net hidden widths; default `c(16L)`),
  `activation` (`"relu"` (default) or `"tanh"`), `epochs` (full-batch
  training epochs, positive integer; default `200L`), `lr` (positive
  Adam learning rate; default `1e-3`), `weight_decay` (non-negative Adam
  L2; default `1e-4`), `min_features` (feature-overlap floor at scoring,
  positive integer; default `3L`), `device` (`"cpu"` (default),
  `"cuda"`, or `"auto"`; `"cuda"`/`"auto"` fall back to CPU with no
  GPU), and `seed` (integer; default `42L`).

## Value

Object of class `coda_deepcoda_model`: a list with `feature_universe`,
`W_bottleneck` (the \\B \times p\\ zero-sum bottleneck weight),
`theta_weights` (exported per-layer `list(W, b)` of the self-explaining
\\\theta\\-net), `bias` (scalar), `activation`, `bottleneck_dim`
(\\B\\), `device` (resolved), `seed`, and `hp`.

## Details

The score is the DeepCoDA logit \\s(x) = \sum_k \theta_k(x)\\ b_k(x) +
c\\ where \\b_k(x) = \sum_j W\_{kj}\\\mathrm{rclr}(x)\_j\\ is the
\\k\\-th zero-sum log-contrast bottleneck (\\\sum_j W\_{kj} = 0\\) and
\\\theta_k(x)\\ the self-explaining per-sample coefficient produced by
the \\\theta\\-net. The zero-sum constraint makes each bottleneck a
scale-invariant log-ratio, so the whole score is exactly per-sample
scale-invariant.

## References

Quinn TP, Nguyen D, Rana S, Gupta S, Venkatesh S. (2020) DeepCoDA:
personalized interpretability for compositional health data. *ICML*,
PMLR 119. arXiv:2006.01392.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
n <- 80; p <- 30; k <- 8
L <- matrix(stats::rnorm(n * p, 4, 0.5), nrow = n,
            dimnames = list(NULL, paste0("miR-", seq_len(p))))
y <- rep(c(0, 1), each = n / 2)
L[y == 1, seq_len(k)] <- L[y == 1, seq_len(k)] + 1.2
X <- exp(L)
model <- fit_coda_deepcoda(X, y)
score_coda_deepcoda(model, X[1, , drop = FALSE])
} # }
```
