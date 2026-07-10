# Fit the SCARF self-supervised contrastive single-sample discriminator

Trains a BatchNorm-free MLP encoder on the per-sample-rCLR,
universe-aligned TRAINING matrix by SELF-SUPERVISED CONTRASTIVE training
(SCARF InfoNCE / NT-Xent, NO labels; via reticulate-python torch, Adam),
FREEZES the encoder, fits a FROZEN linear-probe head on the labeled
training embeddings, EXPORTS the encoder weights and the head \\w, b\\
to R as plain numeric matrices/vectors, and DISCARDS the projection head
and the python module. The head-training embeddings are computed with
the SAME pure-R forward the scorer uses (guaranteeing fit/score
consistency despite float32 training). The fitted model holds no
external pointer and the corruption marginal is NOT carried into it.

## Usage

``` r
fit_ai_scarf(X_train, y_train, meta_train = NULL, hp = list())
```

## Arguments

- X_train:

  Numeric matrix (samples \\\times\\ features) of non-negative
  abundances with unique, non-empty feature names (colnames).

- y_train:

  Numeric / integer 0/1 labels (1 = case), length `nrow(X_train)`, with
  at least one case and one control. Used only to fit the linear-probe
  head; the encoder is trained WITHOUT labels.

- meta_train:

  Optional per-sample metadata. Accepted for interface uniformity and
  otherwise ignored.

- hp:

  Optional list of hyperparameters. Allowed fields: `hidden`
  (positive-integer vector of encoder layer widths up to the embedding;
  the LAST entry is the embedding dim \\d\\; default `c(64L, 32L)`),
  `activation` (`"relu"` (default) or `"tanh"`), `epochs` (number of SSL
  contrastive steps, positive integer; default `200L`), `lr` (positive
  Adam learning rate; default `1e-3`), `weight_decay` (non-negative Adam
  L2 for the SSL stage; default `1e-4`), `corrupt_p` (per-feature
  corruption probability, in \\(0, 1)\\; default `0.6`, the SCARF paper
  value), `temperature` (InfoNCE \\\tau\\, positive; default `1.0`),
  `proj_dim` (projection-head width, positive integer; default `64L`),
  `head_epochs` (linear-probe fixed epochs, positive integer; default
  `200L`), `head_weight_decay` (non-negative Adam L2 for the head;
  default `1e-2`), `min_features` (feature-overlap floor at scoring,
  positive integer; default `3L`), `device` (`"cpu"` (default),
  `"cuda"`, or `"auto"`; `"cuda"`/`"auto"` fall back to CPU with no
  GPU), and `seed` (integer; default `42L`).

## Value

Object of class `ai_scarf_model`: a list with `feature_universe`,
`weights` (exported per-layer `list(W, b)` for the encoder up to the
embedding), `activation`, `head_w` (length-\\d\\ frozen head weight),
`head_b` (scalar frozen head bias), `embedding_dim`, `device`
(resolved), `seed`, and `hp`. The corruption marginal is NOT stored (it
is fit-only).

## Details

The SSL stage learns the encoder by the SCARF InfoNCE loss: each step
builds a random-feature-corrupted view of every training row (corrupted
entries resampled from each feature's empirical training marginal),
embeds the clean row and its corrupted view through encoder + projection
head, and applies the NT-Xent loss (cosine similarity, temperature
\\\tau\\) with the positive pair being the row and its corrupted view
and the negatives being the other batch items. The projection head is
discarded; only the encoder is exported. The head stage fits a single
`Linear(d -> 1)` probe by `head_weight_decay`-regularised BCE on the
frozen encoder embeddings. The score is the linear-probe logit \\s(z) =
w \cdot z + b\\.

## References

Bahri D, Jiang H, Tay Y, Metzler D. (2022) SCARF: Self-Supervised
Contrastive Learning using Random Feature Corruption. *ICLR*.
arXiv:2106.15147.

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
model <- fit_ai_scarf(X, y)
score_ai_scarf(model, X[1, , drop = FALSE])
} # }
```
