# Fit the counterfactual class-conditional VAE single-sample discriminator

Trains a BatchNorm-free CLASS-conditional VAE (encoder \\q(z\|x,y)\\ +
decoder \\p(x\|z,y)\\, both plain MLPs conditioned on a one-hot of the
class label \\y\\) on the per-sample-rCLR, universe-aligned TRAINING
matrix by the reparameterized ELBO (Kingma & Welling 2014; Sohn et al.
2015; via reticulate-python torch, full-batch Adam, minimizing -ELBO),
EXPORTS every weight plus the scalar decoder variance \\\sigma^2\\ to R
as plain numeric matrices/vectors, and DISCARDS the python module. The
fitted model holds no external pointer. The score is the counterfactual
evidence ratio \\\mathrm{ELBO}(x,\mathrm{case}) -
\mathrm{ELBO}(x,\mathrm{control})\\ with the latent MEAN \\z = \mu\\
(deterministic, no RNG), in pure base-R.

## Usage

``` r
fit_cvae(X_train, y_train, meta_train = NULL, hp = list())
```

## Arguments

- X_train:

  Numeric matrix (samples \\\times\\ features) of non-negative
  abundances with unique, non-empty feature names (colnames).

- y_train:

  Numeric / integer 0/1 labels (1 = case), length `nrow(X_train)`, with
  at least one case and one control. Labels condition the VAE at fit
  (the one-hot \\e_y\\).

- meta_train:

  Optional per-sample metadata. Accepted for interface uniformity and
  otherwise ignored.

- hp:

  Optional list of hyperparameters. Allowed fields: `enc_hidden`
  (positive-integer vector of encoder trunk widths; default
  `c(64L, 32L)`), `dec_hidden` (positive-integer vector of decoder trunk
  widths; default `c(32L, 64L)`), `latent_dim` (positive integer latent
  size \\L\\; default `16L`), `activation` (`"relu"` (default) or
  `"tanh"`), `epochs` (number of full-batch training steps, positive
  integer; default `200L`), `lr` (positive Adam learning rate; default
  `1e-3`), `weight_decay` (non-negative Adam L2; default `1e-4`),
  `kl_beta` (non-negative KL coefficient \\\beta\\; default `1.0`),
  `decoder_sigma2` (`NULL` = learned via softplus of a free parameter
  (default), or a positive scalar to fix \\\sigma^2\\), `min_features`
  (feature-overlap floor at scoring, positive integer; default `3L`),
  `device` (`"cpu"` (default), `"cuda"`, or `"auto"`; `"cuda"`/`"auto"`
  fall back to CPU with no GPU), and `seed` (integer; default `42L`).

## Value

Object of class `cvae_model`: a list with `feature_universe`, `weights`
(exported `enc_trunk`/`enc_mu`/`enc_logvar`/ `dec_trunk`/`dec_xhat`),
`sigma2` (scalar decoder variance), `kl_beta`, `latent_dim`,
`activation`, `device` (resolved), `seed`, and `hp`.

## Details

The score is the counterfactual evidence log-likelihood ratio \\s(x) =
\mathrm{ELBO}(x,\mathrm{case}) - \mathrm{ELBO}(x,\mathrm{control})\\.
Departing from the kit-conditional VAE (which conditions on an external
cross-sample covariate), `cvae` conditions on the CLASS label and
evaluates BOTH classes internally, so it needs no kit/cohort input at
score and is single-sample-faithful.

## References

Kingma DP, Welling M. (2014) Auto-Encoding Variational Bayes. *ICLR*.
arXiv:1312.6114.

Sohn K, Lee H, Yan X. (2015) Learning Structured Output Representation
using Deep Conditional Generative Models. *NeurIPS* 28.

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
model <- fit_cvae(X, y)
score_cvae(model, X[1, , drop = FALSE])
} # }
```
