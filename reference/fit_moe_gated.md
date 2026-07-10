# Fit the self-gated mixture-of-experts single-sample discriminator

Trains a LayerNorm MLP encoder, \\E\\ linear logistic expert heads, and
a linear gating network JOINTLY on the per-sample-rCLR, universe-aligned
TRAINING matrix (via reticulate-python torch; Adam, cross-entropy on the
gate-softmax mixture logit), EXPORTS the encoder weights, the LayerNorm
parameters, the \\E\\ expert weight/bias, and the gate weight/bias to R
as plain numeric matrices/vectors, and DISCARDS the python module. The
fitted model holds no external pointer; the score is the pure-R mixture
logit \\s(z) = \sum_e g_e(z)\\h_e(z)\\ with \\z =
\mathrm{enc}(\mathrm{rclr}(x))\\.

The gate is driven by the specimen's OWN embedding (SELF-gated), NOT by
external kit/technology metadata, so the method is
single-sample-deployable (no external covariate at score). See the
method description for this design departure from the manifest's
"kit-gated" name.

## Usage

``` r
fit_moe_gated(X_train, y_train, meta_train = NULL, hp = list())
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
  otherwise ignored (the gate is self-driven, not kit-driven).

- hp:

  Optional list of hyperparameters. Allowed fields: `hidden`
  (positive-integer vector of encoder layer widths; the LAST entry is
  the embedding dim \\d\\; default `c(64L, 32L)`), `activation`
  (`"relu"` (default) or `"tanh"`), `n_experts` (number of experts
  \\E\\, integer \\\geq 2\\; default `4L`), `epochs` (positive integer;
  default `200L`), `lr` (positive Adam learning rate; default `1e-3`),
  `weight_decay` (non-negative Adam L2; default `1e-4`), `min_features`
  (feature-overlap floor at scoring, positive integer; default `3L`),
  `device` (`"cpu"` (default), `"cuda"`, or `"auto"`; `"cuda"`/`"auto"`
  fall back to CPU with no GPU), and `seed` (integer; default `42L`).

## Value

Object of class `moe_gated_model`: a list with `feature_universe`,
`weights` (exported per-layer `list(W, b)` encoder), `layernorm`
(per-layer `list(gamma, beta, eps)`), `activation`, `W_e` (\\E \times
d\\ expert weights), `b_e` (length-\\E\\ expert biases), `W_g` (\\E
\times d\\ gate weights), `b_g` (length-\\E\\ gate biases), `n_experts`,
`embedding_dim`, `device` (resolved), `seed`, and `hp`.

## Details

The score is the convex combination \\s(z) = \sum_e g_e(z)\\h_e(z)\\
where the expert logits are \\h_e(z) = (z W_e^\top + b_e)\_e\\ and the
gate weights are the softmax \\g_e(z) = \mathrm{softmax}(z W_g^\top +
b_g)\_e\\. Because the gate is a softmax, the mixture logit lies within
\\\[\min_e h_e, \max_e h_e\]\\ per row. Larger = more case-like.

## References

Jacobs RA, Jordan MI, Nowlan SJ, Hinton GE. (1991) Adaptive Mixtures of
Local Experts. *Neural Computation* 3(1):79-87. Shazeer N, et al. (2017)
Outrageously Large Neural Networks: The Sparsely-Gated
Mixture-of-Experts Layer. *ICLR*. arXiv:1701.06538.

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
model <- fit_moe_gated(X, y)
score_moe_gated(model, X[1, , drop = FALSE])
} # }
```
