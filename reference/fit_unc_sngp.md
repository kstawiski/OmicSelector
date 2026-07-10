# Fit the spectral-normalized neural Gaussian process discriminator (SNGP)

Trains a spectral-normalized, LayerNorm MLP encoder on the
per-sample-rCLR, universe-aligned TRAINING matrix (via reticulate-python
torch; Adam, cross-entropy on a transient head), EXPORTS the effective
(post-spectral-norm, coefficient-scaled) encoder weights, the LayerNorm
parameters, and a FIXED random-Fourier-feature projection to R as plain
numeric matrices/vectors, DISCARDS the python module and the head, then
fits a closed-form ridge-Laplace logistic posterior (\\\beta\\,
\\\Sigma\\) on the PURE-R RFF features
\\\Phi(\mathrm{enc}(Z\_{\mathrm{train}}))\\ the scorer evaluates
(guaranteeing fit/score consistency). The fitted model holds no external
pointer.

## Usage

``` r
fit_unc_sngp(X_train, y_train, meta_train = NULL, hp = list())
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
  (positive-integer vector of encoder layer widths; the LAST entry is
  the embedding dim \\d\\; default `c(64L, 32L)`), `activation`
  (`"relu"` (default) or `"tanh"`), `epochs` (positive integer; default
  `200L`), `lr` (positive Adam learning rate; default `1e-3`),
  `weight_decay` (non-negative Adam L2; default `1e-4`), `sn_coeff`
  (positive spectral-norm bound \\c\\; default `1.0`), `rff_dim`
  (positive integer number of random Fourier features \\D\\; default
  `128L`), `rff_ls` (positive RBF length-scale \\\ell\\; default `1.0`),
  `ridge` (positive Laplace prior precision \\\lambda\\; default `1.0`),
  `min_features` (feature-overlap floor at scoring, positive integer;
  default `3L`), `device` (`"cpu"` (default), `"cuda"`, or `"auto"`;
  `"cuda"`/`"auto"` fall back to CPU with no GPU), and `seed` (integer;
  default `42L`).

## Value

Object of class `unc_sngp_model`: a list with `feature_universe`,
`weights` (exported per-layer `list(W, b)` effective post-spectral-norm
encoder), `layernorm` (per-layer `list(gamma, beta, eps)`),
`activation`, `W_rff`, `b_rff`, `rff_scale` (\\\sqrt{2/D}\\), `beta`
(Laplace posterior mean), `Sigma` (Laplace posterior covariance),
`Sigma_inv` (precision), `embedding_dim`, `rff_dim`, `device`
(resolved), `seed`, and `hp`.

## Details

The score is the mean-field-adjusted RFF-GP posterior logit \\s(x) =
\Phi(z)\cdot\beta / \sqrt{1 + (\pi/8)\\\Phi(z)^\top\Sigma\\\Phi(z)}\\,
with \\z = \mathrm{enc}(\mathrm{rclr}(x))\\ and \\\Phi(z) =
\sqrt{2/D}\cos(z W\_{\mathrm{rff}}^\top + b\_{\mathrm{rff}})\\. The
denominator is the SNGP distance awareness: the RFF-GP predictive
variance \\\Phi^\top\Sigma\Phi\\ grows with distance from the training
data, shrinking the logit of far-from-data specimens. The Laplace
posterior uses \\\Sigma^{-1} = \lambda I + \sum_i
p_i(1-p_i)\phi_i\phi_i^\top\\.

## References

Liu JZ, Lin Z, Padhy S, Tran D, Bedrax-Weiss T, Lakshminarayanan B.
(2020) Simple and Principled Uncertainty Estimation with Deterministic
Deep Learning via Distance Awareness. *NeurIPS* 33. arXiv:2006.10108.

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
model <- fit_unc_sngp(X, y)
score_unc_sngp(model, X[1, , drop = FALSE])
} # }
```
