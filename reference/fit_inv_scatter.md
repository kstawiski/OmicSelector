# Fit the 1D wavelet-scattering single-sample discriminator

Maps the TRAINING matrix to the per-sample robust CLR over the frozen
feature universe (`colnames(X_train)`), freezes the scattering signal
length to `shape = next_pow2(ncol(X_train))`, computes each training
specimen's float64 1D wavelet scattering coefficients (FORCED
ROW-BY-ROW) via reticulate-python torch + kymatio, freezes a
per-coordinate coefficient standardization, and fits a FROZEN
ridge-logistic head over the standardized coefficients. If `hp$lambda`
is `NULL` the ridge penalty is selected by TRAINING-ONLY deterministic
class-stratified k-fold CV (no RNG). The fitted model holds NO external
pointer – score rebuilds `ScatteringTorch1D` from the stored config.

## Usage

``` r
fit_inv_scatter(X_train, y_train, meta_train = NULL, hp = list())
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

  Optional list of hyperparameters. Allowed fields: `J` (scattering
  scale, positive integer; default `3L`), `Q` (wavelets per octave,
  positive integer; default `4L`), `max_order` (scattering order, `1L`
  or `2L`; default `2L`), `lambda` (ridge penalty for the logistic head;
  `NULL` (default) selects it by training-only deterministic k-fold CV,
  otherwise a non-negative finite number), `min_features`
  (feature-overlap floor at scoring, positive integer; default `3L`),
  `device` (`"cpu"` (default), `"cuda"`, or `"auto"`; `"cuda"`/`"auto"`
  fall back to CPU with no GPU), and `seed` (integer; default `42L`).

## Value

Object of class `inv_scatter_model`: a list with `feature_universe`,
`J`, `Q`, `shape`, `max_order`, `coef_mean`, `coef_sd` (frozen
standardization), `intercept`, `weights` (frozen linear head), `lambda`
(resolved ridge penalty), `device` (resolved), `seed`, and `hp`. No
python pointer.

## References

Anden J, Mallat S. (2014) Deep Scattering Spectrum. *IEEE TSP*
62(16):4114-4128.

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
model <- fit_inv_scatter(X, y)
score_inv_scatter(model, X[1, , drop = FALSE])
} # }
```
