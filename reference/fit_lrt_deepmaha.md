# Fit the deep-feature class-conditional Mahalanobis LRT discriminator

Trains a BatchNorm-free MLP embedding on the per-sample-rCLR,
universe-aligned TRAINING matrix (via reticulate-python torch; Adam,
cross-entropy), EXPORTS the weights up to the penultimate embedding
layer to R as plain numeric matrices/vectors, DISCARDS the python module
and the 2-logit head, then estimates the frozen per-class Gaussian means
and a single TIED ridge-regularised covariance in the embedding space
using the PURE-R forward the scorer evaluates (guaranteeing fit/score
consistency despite float32 training). The fitted model holds no
external pointer.

## Usage

``` r
fit_lrt_deepmaha(X_train, y_train, meta_train = NULL, hp = list())
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

  Optional list of hyperparameters. Allowed fields: `hidden`
  (positive-integer vector of layer widths up to the embedding; the LAST
  entry is the embedding dim \\d\\; default `c(64L, 32L)`), `activation`
  (`"relu"` (default) or `"tanh"`), `epochs` (positive integer; default
  `200L`), `lr` (positive Adam learning rate; default `1e-3`),
  `weight_decay` (non-negative Adam L2; default `1e-4`), `batch_size`
  (positive integer minibatch; default `32L`), `ridge` (positive
  covariance ridge factor; default `1e-3`), `min_features`
  (feature-overlap floor at scoring, positive integer; default `3L`),
  `device` (`"cpu"` (default), `"cuda"`, or `"auto"`; `"cuda"`/`"auto"`
  fall back to CPU with no GPU), and `seed` (integer; default `42L`).

## Value

Object of class `lrt_deepmaha_model`: a list with `feature_universe`,
`weights` (exported per-layer `list(W, b)` up to the embedding),
`activation`, `mu_case`, `mu_ctrl`, `Sigma_reg_inv` (inverse of the
ridge-regularised tied covariance), `embedding_dim`, `device`
(resolved), `seed`, and `hp`.

## Details

The score is the class-conditional Gaussian log-likelihood ratio with a
shared covariance: \\s(z) = -\tfrac12 (z-\mu\_{\mathrm{case}})^\top
\Sigma^{-1} (z-\mu\_{\mathrm{case}}) + \tfrac12
(z-\mu\_{\mathrm{ctrl}})^\top \Sigma^{-1} (z-\mu\_{\mathrm{ctrl}})\\,
exactly linear in the embedding \\z\\ (LDA in embedding space). The tied
covariance is the pooled within-class estimator \\\Sigma = ((n_1-1)S_1 +
(n_0-1)S_0)/(n_1+n_0-2)\\, ridge-regularised as \\\Sigma +
\lambda\\\overline{\mathrm{diag}(\Sigma)}\\I\\ with \\\lambda = \\
`hp$ridge`, with a positive-definite floor if needed. This uses the
unbiased pooled denominator; Lee 2018 displays the \\1/N\\ MLE form,
which differs only by a constant scaling of the tied covariance and so
leaves the per-specimen score ranking (and AUC) unchanged.

## References

Lee K, Lee K, Lee H, Shin J. (2018) A Simple Unified Framework for
Detecting Out-of-Distribution Samples and Adversarial Attacks. *NeurIPS*
31. arXiv:1807.03888.

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
model <- fit_lrt_deepmaha(X, y)
score_lrt_deepmaha(model, X[1, , drop = FALSE])
} # }
```
