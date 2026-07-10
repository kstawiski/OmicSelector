# Fit the V-REx cross-cohort transfer single-sample discriminator

Trains a BatchNorm-free MLP encoder + a linear head \\\mathrm{Linear}(d
\to 1)\\ JOINTLY across training environments (cohorts) by the V-REx
objective (Krueger 2021; via reticulate-python torch, Adam, mean
per-environment risk plus a population-variance-of-risks penalty),
EXPORTS the encoder weights up to and including the embedding layer AND
the head \\(w, b)\\ to R as plain numeric arrays, and DISCARDS the
python module. The fitted model holds no external pointer; scoring is
pure base-R.

Environments are the distinct non-missing labels in
`meta_train[[cohort_col]]` (`cohort_col` auto-detected when `NULL`).
Only environments carrying both a case and a control are retained; if
fewer than two valid environments remain (including `NULL`/absent
`meta_train`), fitting falls back to plain ERM (a single pooled risk,
`vrex_beta` inert) and `model$n_environments` is set to `1`.

## Usage

``` r
fit_dro_vrex(X_train, y_train, meta_train = NULL, hp = list())
```

## Arguments

- X_train:

  Numeric matrix (samples \\\times\\ features) of non-negative
  abundances with unique, non-empty feature names (colnames).

- y_train:

  Numeric / integer 0/1 labels (1 = case), length `nrow(X_train)`, with
  at least one case and one control.

- meta_train:

  Optional per-sample metadata carrying a cohort/environment column.
  When it contains the resolved `cohort_col` with at least two
  both-class cohorts, V-REx trains across those environments; otherwise
  the fit degenerates to ERM. Must have one row per row of `X_train`.

- hp:

  Optional list of hyperparameters. Allowed fields: `hidden`
  (positive-integer vector of encoder layer widths up to the embedding;
  the LAST entry is the embedding dim \\d\\; default `c(64L, 32L)`),
  `activation` (`"relu"` (default) or `"tanh"`), `epochs` (training
  steps, positive integer; default `200L`), `lr` (positive Adam learning
  rate; default `1e-3`), `weight_decay` (non-negative Adam L2; default
  `1e-4`), `vrex_beta` (non-negative V-REx variance-penalty coefficient;
  default `1.0`), `cohort_col` (`NULL` (default, auto-detect) or a
  single character string naming the environment column in
  `meta_train`), `min_features` (feature-overlap floor at scoring,
  positive integer; default `3L`), `device` (`"cpu"` (default),
  `"cuda"`, or `"auto"`), and `seed` (integer; default `42L`).

## Value

Object of class `dro_vrex_model`: a list with `feature_universe`,
`weights` (exported per-layer `list(W, b)` up to the embedding),
`activation`, `head_w` (length-\\d\\ frozen head weight), `head_b`
(scalar frozen head bias), `embedding_dim`, `n_environments` (number of
valid both-class environments engaged; `1` means ERM fallback), `device`
(resolved), `seed`, and `hp`.

## References

Krueger D, et al. (2021) Out-of-Distribution Generalization via Risk
Extrapolation (REx). *ICML*. arXiv:2003.00688.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
n <- 120; p <- 30; k <- 8
L <- matrix(stats::rnorm(n * p, 4, 0.5), nrow = n,
            dimnames = list(NULL, paste0("miR-", seq_len(p))))
y <- rep(c(0, 1), each = n / 2)
L[y == 1, seq_len(k)] <- L[y == 1, seq_len(k)] + 1.2
X <- exp(L)
meta <- data.frame(accession = rep(paste0("GSE", 1:3), length.out = n))
model <- fit_dro_vrex(X, y, meta_train = meta)
score_dro_vrex(model, X[1, , drop = FALSE])
} # }
```
