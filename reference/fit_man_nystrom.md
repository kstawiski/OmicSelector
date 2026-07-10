# Fit a Nystrom diffusion-map discriminator on frozen landmarks

Learns a single-sample circulating-miRNA discriminator by selecting
frozen training landmarks, building a diffusion-map geometry on their
per-sample rCLR coordinates, embedding every training specimen by the
Nystrom out-of-sample formula, and fitting a frozen linear head. The
rCLR transform is computed independently for each specimen: positive
entries are log transformed and centered by that specimen's own mean
log-positive abundance, while zero entries remain at 0.

## Usage

``` r
fit_man_nystrom(X_train, y_train, meta_train = NULL, hp = list())
```

## Arguments

- X_train:

  Numeric matrix (samples \\\times\\ features) of non-negative abundance
  values. Columns must be uniquely named feature ids.

- y_train:

  Integer/numeric 0/1 labels aligned to `X_train`; 1 is case/disease and
  0 is control.

- meta_train:

  Optional per-sample metadata. Accepted for interface uniformity and
  ignored by this method.

- hp:

  List of frozen hyperparameters: `n_landmarks` positive integer
  landmark cap (default `min(200L, nrow(X_train))`); `n_components`
  positive integer diffusion-coordinate cap (default `10L`); `t`
  positive integer diffusion time (default `1L`); `gamma` positive RBF
  bandwidth, or `NULL` for the landmark median heuristic (default
  `NULL`); `eig_floor` non-negative eigenvalue floor for retained
  non-trivial modes (default `1e-6`); `min_features` positive integer
  scoring overlap floor (default `3L`); `eps` positive numerical floor
  (default `1e-8`); and `seed` non-negative integer used only for
  deterministic landmark subsampling while restoring the global RNG
  state (default `1L`).

## Value

A plain list of class `man_nystrom_model` containing the frozen training
`feature_universe`, rCLR `landmarks`, `landmark_idx`, frozen `gamma`,
diffusion settings, cached full-universe `geometry`, frozen linear
`head`, and resolved `hp`.

## Details

The landmark graph is learned once from training data. For landmarks
\\L_j\\, the RBF kernel \\W\_{ij}=\exp(-\gamma\\L_i-L_j\\^2)\\ is
symmetrically degree-normalized and eigendecomposed. The trivial top
diffusion pair is dropped; retained right eigenvectors of the
row-stochastic transition matrix define the diffusion coordinates. A new
specimen \\z\\ is embedded from only its own affinities to the frozen
landmarks: \$\$\Psi_k(z)=\lambda_k^{t-1}\sum_j p_j(z)\psi\_{jk},\$\$
where \\p_j(z)\\ is the row-normalized RBF affinity of \\z\\ to landmark
\\j\\. These coordinates are standardized by frozen training
means/scales and passed to a frozen binomial-GLM head, with a
diagonal-LDA mean-difference fallback if the GLM fit is numerically
unsafe. Larger scores are oriented to be more case-like.

No scored-batch statistic is used. At scoring, the present feature set
is `intersect(model$feature_universe, colnames(X))`; the landmark
geometry is re-derived from frozen landmark coordinates restricted to
those features, and each scored row contributes only its own per-sample
rCLR vector and its own affinity row to the frozen landmarks. If fewer
than `hp$min_features` model features are present, the documented
neutral score `0` is returned for every row.

## References

Coifman RR, Lafon S. (2006) Diffusion maps. *Applied and Computational
Harmonic Analysis* 21: 5-30.

Bengio Y, Paiement JF, Vincent P, Delalleau O, Le Roux N, Ouimet M.
(2004) Out-of-sample extensions for LLE, Isomap, MDS, Eigenmaps, and
Spectral Clustering. *Advances in Neural Information Processing Systems*
16.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
X <- matrix(stats::rgamma(60 * 30, shape = 3), nrow = 60,
            dimnames = list(NULL, paste0("miR-", seq_len(30))))
y <- rep(c(0, 1), each = 30)
X[y == 1, 1:6] <- X[y == 1, 1:6] + 20
X[y == 0, 7:12] <- X[y == 0, 7:12] + 20
model <- fit_man_nystrom(X, y)
score <- score_man_nystrom(model, X)
} # }
```
