# Fit the GASF-image CNN single-sample discriminator

Maps the TRAINING matrix to the per-sample robust CLR over the frozen
feature universe (`colnames(X_train)`), computes the FROZEN per-feature
GASF bounds, renders each training specimen to a GASF image, and trains
a small Conv-BN-ReLU-MaxPool CNN end-to-end by full-batch
`BCEWithLogitsLoss` (Adam, no shuffle, seed-before-build) via
reticulate-python torch + torchvision. After training the net is put in
[`eval()`](https://rdrr.io/r/base/eval.html) mode (frozen BatchNorm) and
its FULL `state_dict` is EXPORTED to R as named float64 numeric arrays
(incl. BN running stats); the python module is DISCARDED. The fitted
model holds NO external pointer – score rebuilds the CNN from the R
state.

## Usage

``` r
fit_img_gasfcnn(X_train, y_train, meta_train = NULL, hp = list())
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

  Optional list of hyperparameters. Allowed fields: `channels`
  (positive-integer vector of per-block conv channel widths; default
  `c(8L, 16L)` – two blocks), `epochs` (full-batch training epochs,
  positive integer; default `200L`), `lr` (positive Adam learning rate;
  default `1e-3`), `weight_decay` (non-negative Adam L2; default
  `1e-4`), `min_features` (feature-overlap floor at scoring, positive
  integer; default `3L`), `device` (`"cpu"` (default), `"cuda"`, or
  `"auto"`; `"cuda"`/`"auto"` fall back to CPU with no GPU), `seed`
  (integer; default `42L`), and `image_size` (fixed GASF image side
  length, integer \\\ge 2\\; default `128L`). The frozen-order rCLR
  profile is PAA-reduced to length \\L = \min(p, image\\size)\\ before
  imaging, so the GASF image is at most `image_size` x `image_size`
  regardless of the feature count \\p\\; when \\p \le image\\size\\ no
  reduction is applied.

## Value

Object of class `img_gasfcnn_model`: a list with `feature_universe`,
`state_dict` (named float64 arrays + BN counters), `channels`,
`gasf_lo`, `gasf_hi` (frozen per-coordinate GASF bounds, length
`paa_len`), `image_size` (the PAA target side length) and `paa_len`
(`= min(p, image_size)`, the actual GASF image side), `device`
(resolved), `seed`, and `hp`. No python pointer.

## References

Wang Z, Oates T. (2015) Imaging Time-Series to Improve Classification
and Imputation. *IJCAI*; arXiv:1506.00327.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
n <- 80; p <- 16; k <- 6
L <- matrix(stats::rnorm(n * p, 4, 0.5), nrow = n,
            dimnames = list(NULL, paste0("miR-", seq_len(p))))
y <- rep(c(0, 1), each = n / 2)
L[y == 1, seq_len(k)] <- L[y == 1, seq_len(k)] + 1.2
X <- exp(L)
model <- fit_img_gasfcnn(X, y)
score_img_gasfcnn(model, X[1, , drop = FALSE])
} # }
```
