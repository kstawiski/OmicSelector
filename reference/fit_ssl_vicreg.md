# Fit the VICReg self-supervised single-sample discriminator

Trains a BatchNorm-free MLP encoder on the per-sample-rCLR,
universe-aligned TRAINING matrix by SELF-SUPERVISED VICReg
regularization (Bardes 2022; via reticulate-python torch, Adam,
variance-invariance-covariance loss on two augmented views – NO labels),
EXPORTS the encoder weights up to and including the embedding layer to R
as plain numeric matrices/vectors, DISCARDS the expander and the python
module, then fits a FROZEN linear-probe head (\\w, b\\) on the labeled
training embeddings computed by the PURE-R forward the scorer uses
(guaranteeing fit/score consistency despite float32 training). The
fitted model holds no external pointer.

## Usage

``` r
fit_ssl_vicreg(X_train, y_train, meta_train = NULL, hp = list())
```

## Arguments

- X_train:

  Numeric matrix (samples \\\times\\ features) of non-negative
  abundances with unique, non-empty feature names (colnames).

- y_train:

  Numeric / integer 0/1 labels (1 = case), length `nrow(X_train)`, with
  at least one case and one control. Labels are used ONLY for the frozen
  linear-probe head; the encoder is trained without labels.

- meta_train:

  Optional per-sample metadata. Accepted for interface uniformity and
  otherwise ignored.

- hp:

  Optional list of hyperparameters. Allowed fields: `hidden`
  (positive-integer vector of encoder layer widths up to the embedding;
  the LAST entry is the embedding dim \\d\\; default `c(64L, 32L)`),
  `activation` (`"relu"` (default) or `"tanh"`), `epochs` (number of SSL
  training steps, positive integer; default `200L`), `lr` (positive Adam
  learning rate; default `1e-3`), `weight_decay` (non-negative
  encoder/expander Adam L2; default `1e-4`), `mask_p` (per-view
  feature-mask fraction in (0, 1); default `0.3`), `noise_sd`
  (non-negative rCLR Gaussian-jitter sd; default `0.1`), `expander_dim`
  (positive-integer expander/projector width; default `64L`),
  `head_epochs` (linear-probe fixed epochs, positive integer; default
  `200L`), `head_weight_decay` (non-negative head Adam L2; default
  `1e-2`), `min_features` (feature-overlap floor at scoring, positive
  integer; default `3L`), `device` (`"cpu"` (default), `"cuda"`, or
  `"auto"`; `"cuda"`/`"auto"` fall back to CPU with no GPU), and `seed`
  (integer; default `42L`).

## Value

Object of class `ssl_vicreg_model`: a list with `feature_universe`,
`weights` (exported per-layer `list(W, b)` up to the embedding),
`activation`, `head_w` (length-\\d\\ frozen head weight), `head_b`
(scalar frozen head bias), `embedding_dim`, `device` (resolved), `seed`,
and `hp`.

## Details

The score is the frozen linear-probe logit \\s(z) = w^\top z + b\\ on
the encoder embedding \\z\\. The distinction from the cross-entropy
`lrt-deepmaha` and prototypical `proto-net` siblings is the TRAINING:
the encoder is learned with NO labels by the VICReg loss (Bardes 2022
eq. 6, \\\lambda = 25, \mu = 25, \nu = 1\\) on two augmented views, and
only a small linear head reads the labels.

## References

Bardes A, Ponce J, LeCun Y. (2022) VICReg:
Variance-Invariance-Covariance Regularization for Self-Supervised
Learning. *ICLR*. arXiv:2105.04906.

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
model <- fit_ssl_vicreg(X, y)
score_ssl_vicreg(model, X[1, , drop = FALSE])
} # }
```
