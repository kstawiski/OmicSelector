# Fit the prototypical-network single-sample discriminator

Trains a BatchNorm-free MLP embedding on the per-sample-rCLR,
universe-aligned TRAINING matrix by EPISODIC PROTOTYPICAL training
(Snell 2017; via reticulate-python torch, Adam, prototypical loss – NOT
cross-entropy), EXPORTS the weights up to and including the embedding
layer to R as plain numeric matrices/vectors, DISCARDS the python
module, then computes the frozen per-class PROTOTYPES (class-mean
embeddings) using the PURE-R forward over ALL training rows
(guaranteeing fit/score consistency despite float32 training). The
fitted model holds no external pointer.

## Usage

``` r
fit_proto_net(X_train, y_train, meta_train = NULL, hp = list())
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
  (`"relu"` (default) or `"tanh"`), `epochs` (number of training
  episodes, positive integer; default `200L`), `lr` (positive Adam
  learning rate; default `1e-3`), `weight_decay` (non-negative Adam L2;
  default `1e-4`), `n_support` (per-class support size per episode,
  positive integer; default `5L`), `n_query` (per-class query size per
  episode, positive integer; default `5L`), `min_features`
  (feature-overlap floor at scoring, positive integer; default `3L`),
  `device` (`"cpu"` (default), `"cuda"`, or `"auto"`; `"cuda"`/`"auto"`
  fall back to CPU with no GPU), and `seed` (integer; default `42L`).
  `n_support`/`n_query` are clamped per class so each episode is
  feasible (disjoint support/query, \\\geq 1\\ of each per class).

## Value

Object of class `proto_net_model`: a list with `feature_universe`,
`weights` (exported per-layer `list(W, b)` up to the embedding),
`activation`, `proto_case`, `proto_ctrl` (length-\\d\\ frozen
prototypes), `embedding_dim`, `device` (resolved), `seed`, and `hp`.

## Details

The score is the squared-Euclidean log-likelihood ratio to the
prototypes, \\s(z) = -\lVert z - \mu\_{\mathrm{case}}\rVert^2 + \lVert
z - \mu\_{\mathrm{ctrl}}\rVert^2\\, exactly the deepmaha score with
\\\Sigma = I\\. The distinction is the embedding: it is trained by the
prototypical loss (CE over the softmax of negative squared-Euclidean
distances to episodic class-mean prototypes), not a cross-entropy
classification head.

## References

Snell J, Swersky K, Zemel R. (2017) Prototypical Networks for Few-shot
Learning. *NeurIPS* 30. arXiv:1703.05175.

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
model <- fit_proto_net(X, y)
score_proto_net(model, X[1, , drop = FALSE])
} # }
```
