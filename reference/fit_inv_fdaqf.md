# Fit quantile-function FDA (inv-fdaqf) within-sample discriminator

Learns a frozen, single-sample discriminator from a functional-PCA
(FPCA) decomposition of a specimen's compositional *quantile function*
(QF). Each specimen is reduced to a single monotone curve – the sorted
robust-CLR (rCLR) abundance profile, resampled to a fixed length – and
that curve is projected onto a frozen orthonormal FPCA eigenbasis, whose
coordinates feed a frozen ridge-regularised linear-discriminant (LDA)
head. The construction has four stages, all estimated from training data
only (the QF grid is a fixed data-free constant):

1.  **Single-specimen quantile function.** A specimen's non-negative
    abundances are mapped to a within-sample robust-CLR vector: over the
    specimen's own positive support \\P=\\j: v_j\>0\\\\, \\z_j = \log
    v_j - \mathrm{mean}\_{k\in P}\log v_k\\. Structural zeros
    (\\v_j=0\\) are excluded because they encode detection-rate /
    sequencing depth, a batch artifact rather than shape (only
    structural zeros, keyed on \\v_j=0\\, are dropped; a present
    coordinate with \\z_j=0\\, e.g. a tied / constant count profile, is
    kept). The present rCLR coordinates \\z\_{pos}\\ are resampled to a
    fixed length \\L\\ by evaluating their empirical quantile function
    at \\u\_\ell=(\ell-0.5)/L\\, \\\ell=1,\dots,L\\
    ([`stats::quantile`](https://rdrr.io/r/stats/quantile.html),
    `type = 7`). The result \\q\in\mathbb{R}^L\\ is a monotone
    non-decreasing curve of fixed length \\L\\ regardless of the
    specimen's support size, and *exactly* invariant to per-sample
    positive scaling (rCLR removes the per-sample geometric mean, so
    multiplying \\v\\ by any \\c\>0\\ leaves \\z\_{pos}\\, hence \\q\\,
    unchanged). This QF substrate is shared with the
    curvature-scale-space (inv-css) and optimal-transport profile
    methods, but the mechanism here is a distinct functional-PCA
    projection, not a curvature descriptor and not an optimal-transport
    embedding.

2.  **Frozen functional-PCA basis (training only).** The training QF
    matrix \\G\\ (\\n\_{kept}\times L\\, one row per retained specimen
    QF) is centred by its functional mean \\\bar
    q=\mathrm{colMeans}(G)\\: \\G_c = G - \mathbf{1}\bar q^\top\\. The
    top-\\K\\ right-singular vectors of \\G_c\\ (base-R `svd`, no fda
    dependency) form the frozen orthonormal loading basis \\V\\
    (\\L\times K\\, with \\V^\top V = I_K\\ to machine precision).
    \\K=\\`n_components`, capped at the numerical rank of \\G_c\\: if
    fewer than \\K\\ positive singular values exist only those are kept
    and the realised \\K\\ is recorded (`model$fpca$K`). SVD column
    signs (and, for repeated singular values, the in-subspace rotation)
    are arbitrary but *frozen* in the model and absorbed by the
    downstream linear head, so the final scores are sign-robust; the
    frozen \\V\\ orientation is nonetheless LAPACK-build dependent (a
    reproducibility caveat shared with other frozen eigenbasis methods
    such as gsp-gft).

3.  **Projection scores.** A specimen's FPCA score vector is the
    projection of its centred QF onto the frozen basis, \\s = V^\top
    (q - \bar q)\in\mathbb{R}^K\\. It uses only the specimen's own QF
    and the frozen \\\bar q\\, \\V\\.

4.  **Frozen standardize + ridge-LDA head.** The \\K\\-dimensional score
    vectors are standardised by a frozen `center`/`scale` (each
    component's training mean and standard deviation, the latter floored
    at `eps`; exactly-zero-variance components are deactivated and
    contribute `0`). A ridge-regularised LDA head is then fit: pooled
    within-class covariance \\S_w\\, ridged to positive definiteness
    (\\S_w + (\mathrm{shrink}\cdot\overline{\mathrm{diag}} +
    \mathrm{eps})I\\, the ridge raised until the Cholesky succeeds), \\w
    = S_w^{-1}(\mu\_{case}-\mu\_{ctrl})\\ and \\b = -\tfrac12
    (\mu\_{case}+\mu\_{ctrl})^\top w\\. The head orients *larger = more
    case-like*.

Disease reshapes the *distributional shape* of the compositional profile
(e.g. a heavier upper tail in the rCLR values); the FPCA score vector is
a compact, batch-robust, scale-free summary of that shape along the
dominant training modes of variation. No test data and no test-time
batch statistic are used: the rCLR is strictly within-sample, and the QF
grid, functional mean, eigenbasis, standardisation and head are all
frozen here from training only.

## Usage

``` r
fit_inv_fdaqf(X_train, y_train, meta_train = NULL, hp = list())
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

  List of frozen hyperparameters (exact-name reads, strict
  allowed-list): `resample_len` fixed QF length \\L\\ (default `64L`;
  integer \\\ge 16\\), `n_components` number of frozen FPCA components
  \\K\\ (default `6L`; integer \\\ge 1\\; the realised \\K\\ is reduced
  at fit time when the training centred QFs span a subspace of rank \\\<
  K\\), `shrink` ridge fraction of the mean covariance diagonal (default
  `0.1`; \\\ge 0\\), `min_features` per-specimen floor on the number of
  present rCLR coordinates (default `8L`; integer \\\ge 3\\), `eps`
  positive standard-deviation / ridge floor (default `1e-6`), and `seed`
  (default `1L`; accepted for interface uniformity only – the method
  consumes no randomness).

## Value

A plain list of class `inv_fdaqf_model` containing `feature_universe`,
the frozen `fpca` basis (functional mean `qbar`, orthonormal loadings
`V`, realised `K`, numerical `rank`, training `singular_values`),
`descriptor_dim` (\\= K\\), the frozen `head` (`center`, `scale`,
`active`, `w`, `b`, `ridge`), and the resolved `hp`.

## Details

The per-specimen score (see
[`score_inv_fdaqf`](https://kstawiski.github.io/OmicSelector/reference/score_inv_fdaqf.md))
is the frozen linear predictor on the standardised FPCA score vector,
\$\$S(x) = b + \sum\_{c} w_c \\ \frac{s_c(x) - \mathrm{center}\_c}
{\mathrm{scale}\_c}, \qquad s(x) = V^\top (q(x) - \bar q),\$\$ where
\\q(x)\\ is the specimen's quantile function. Every quantity except
\\q(x)\\ is frozen at fit time, and \\q(x)\\ depends only on \\x\\'s own
rCLR profile, so the score is computable from a single specimen and is
exactly invariant to per-sample positive scaling. Although it reads the
same sorted rCLR substrate as the curvature-scale-space and
optimal-transport profile methods, the representation here is a distinct
functional-PCA projection.

## References

Ramsay JO, Silverman BW. (2005) Functional Data Analysis, 2nd ed.
Springer. (Functional principal-component analysis.)

Bolstad BM, Irizarry RA, Astrand M, Speed TP. (2003) A comparison of
normalization methods for high density oligonucleotide array data based
on variance and bias. *Bioinformatics* 19(2): 185-193. (Quantile
representation of a sample distribution.)

Martino C, Morton JT, Marotz CA, et al. (2019) A novel sparse
compositional technique reveals microbial perturbations. *mSystems*
4(1): e00016-19. (Robust CLR.)

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
p <- 60
features <- paste0("miR-", sprintf("%03d", seq_len(p)))
n <- 160
y <- rep(c(0, 1), each = n / 2)
X <- matrix(0, n, p, dimnames = list(NULL, features))
rank_frac <- (seq_len(p) - 1) / (p - 1)
ctrl_mu <- 2 + 3 * rank_frac                               # near-linear QF
case_mu <- 2 + 3 * rank_frac^3                             # heavier upper tail
for (i in seq_len(n)) {
  mu <- if (y[i] == 1) case_mu else ctrl_mu
  X[i, ] <- exp(stats::rnorm(p, mean = mu, sd = 0.2))
}
model <- fit_inv_fdaqf(X, y)
score <- score_inv_fdaqf(model, X)
} # }
```
