# Fit a graph-Fourier discriminator on a frozen co-expression graph

Learns a single-sample circulating-miRNA discriminator by building a
frozen feature-feature co-expression graph from the training set,
extracting the low-frequency graph-Fourier modes of its
symmetric-normalized Laplacian, and fitting a frozen linear head on
per-specimen graph-spectral coordinates. Each specimen is represented by
a robust centered log-ratio (rCLR) transform computed from that specimen
alone: positive entries are log-transformed and centered by the mean log
abundance of the specimen's own positive entries, while zeros remain
neutral at 0.

## Usage

``` r
fit_gsp_gft(X_train, y_train, meta_train = NULL, hp = list())
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

  List of frozen hyperparameters: `cor_method`, one of `"spearman"` or
  `"pearson"` (default `"spearman"`); `graph_k`, positive integer number
  of nearest neighbours retained per feature before symmetrization
  (default `min(10L, D - 1L)`); `n_modes`, positive integer cap on
  retained low-frequency graph Fourier modes (default `min(10L, D)`);
  `min_features`, positive integer scoring overlap floor (default `3L`);
  and `eps`, positive finite numerical floor (default `1e-8`).

## Value

A plain list of class `gsp_gft_model` containing the frozen
`feature_universe`, low-frequency basis `U_low`, retained eigenvalues
`eigvals_low`, graph settings, frozen linear head, and the resolved
hyperparameters.

## Details

The co-expression graph is learned from training rCLR signals only. The
absolute training feature-feature correlation matrix is kNN-sparsified
with a deterministic feature-index tie break, symmetrized by union, and
converted to the symmetric-normalized graph Laplacian \$\$L = I -
D^{-1/2} W D^{-1/2}.\$\$ Its eigenvectors, ordered by increasing
eigenvalue, define graph-Fourier coordinates. The retained low-frequency
coordinates are standardized by frozen training means and scales and
then passed to a frozen linear head. A binomial GLM is used when it is
finite and converged; otherwise the fit falls back to a diagonal-LDA
mean-difference head in the same standardized coordinates. Larger scores
are oriented to be more case-like.

At inference, no statistic is estimated from the scoring matrix. The
graph, basis, coordinate center/scale, and head coefficients are all
frozen at fit time from the training data. Each row's rCLR vector uses
only that row's own positive present features, with absent and zero
features encoded as 0.

## References

Shuman DI, Narang SK, Frossard P, Ortega A, Vandergheynst P. (2013) The
emerging field of signal processing on graphs. *IEEE Signal Processing
Magazine* 30: 83-98.

Chung FRK. (1997) *Spectral Graph Theory*. American Mathematical
Society.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
X <- matrix(stats::rgamma(60 * 30, shape = 3), nrow = 60,
            dimnames = list(NULL, paste0("miR-", seq_len(30))))
y <- rep(c(0, 1), each = 30)
X[y == 1, 1:6] <- X[y == 1, 1:6] + 20
model <- fit_gsp_gft(X, y)
score <- score_gsp_gft(model, X)
} # }
```
