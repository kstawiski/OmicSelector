# Fit curvature-scale-space (inv-css) within-sample discriminator

Learns a frozen, single-sample discriminator from the *multi-scale
differential curvature* of a specimen's compositional shape. Each
specimen is reduced to a single 1-D curve – the sorted robust-CLR (rCLR)
abundance profile – and that curve is summarised by a
Mokhtarian-Mackworth curvature-scale-space (CSS) descriptor, which is
then fed to a frozen linear (LDA) head. The construction has three
stages, all estimated from training data only (or, like the scale ladder
and kernels, fixed data-free constants):

1.  **Single-specimen 1-D curve.** A specimen's non-negative abundances
    are mapped to a within-sample robust-CLR vector: over the specimen's
    own positive support \\P=\\j: v_j\>0\\\\, \\z_j = \log v_j -
    \mathrm{mean}\_{k\in P}\log v_k\\; structural zeros (\\v_j=0\\) are
    excluded because they encode detection-rate / sequencing depth, a
    batch artifact rather than shape (only structural zeros, keyed on
    \\v_j=0\\, are dropped; a present coordinate with \\z_j=0\\ is
    kept). The present rCLR coordinates \\z\_{pos}\\ are resampled to a
    fixed length \\L\\ by evaluating their empirical quantile function
    at \\u\_\ell=(\ell-0.5)/L\\, \\\ell=1,\dots,L\\
    ([`stats::quantile`](https://rdrr.io/r/stats/quantile.html),
    `type = 7`). The result \\g\in\mathbb{R}^L\\ is a monotone
    non-decreasing curve, approximately mean-zero (the midpoint quantile
    resample of the exactly-mean-zero \\z\_{pos}\\ need not preserve the
    mean exactly; the curvature descriptors are in any case invariant to
    a vertical shift of \\g\\), and *exactly* invariant to per-sample
    positive scaling (rCLR removes the per-sample geometric mean, so
    multiplying \\v\\ by any \\c\>0\\ leaves \\z\_{pos}\\, hence \\g\\,
    unchanged).

2.  **Curvature-scale-space descriptor.** A frozen geometric ladder of
    \\S\\ scales \\\sigma_s\\ runs from \\\sigma\_{min}\\ to
    \\\sigma\_{max}\\. At each scale \\g\\ is convolved with a discrete,
    sum-1 Gaussian kernel (half-width \\h=\lceil 3\sigma_s\rceil\\,
    symmetric/"reflect" boundary padding) to give \\g_s\\; unit-grid
    central and second differences \\g'\_s, g''\_s\\ give the planar
    curve curvature \\\kappa_s = g''\_s / (1 + g'^2_s)^{3/2}\\ at the
    interior grid points, summarised over the **padding-free
    sub-interval** only (indices whose smoothing window lies entirely
    within the real signal, excluding the reflect-padding boundary
    artefact that a monotone curve's mirrored slope would otherwise
    inject). Each scale contributes four descriptors – total absolute
    curvature \\\sum\|\kappa_s\|\\, curvature energy \\\sum\kappa_s^2\\,
    peak curvature \\\max\|\kappa_s\|\\, and the CSS inflection count
    (number of sign changes of \\\kappa_s\\) – concatenated over scales
    into \\\varphi\in\mathbb{R}^{4S}\\. The scale ladder, kernels and
    difference stencils are data-free, so the descriptor is defined
    identically at fit and at scoring.

3.  **Frozen standardize + LDA head.** Training descriptors are
    standardised by a frozen `center`/`scale` (each descriptor's
    training mean and standard deviation, the latter floored at `eps`;
    exactly-zero-variance descriptors are deactivated and contribute
    `0`). A ridge-regularised linear-discriminant head is fit on the
    standardised descriptors: pooled within-class covariance \\S_w\\,
    ridged to positive definiteness (\\S_w +
    (\mathrm{shrink}\cdot\overline{\mathrm{diag}} + \mathrm{eps})I\\,
    the ridge raised until the Cholesky succeeds), \\w =
    S_w^{-1}(\mu\_{case}-\mu\_{ctrl})\\ and \\b = -\tfrac12
    (\mu\_{case}+\mu\_{ctrl})^\top w\\. The head orients *larger = more
    case-like*.

Disease reshapes the *shape* of the sorted compositional profile (a
sharp "shoulder" or kink versus a smooth ramp); the CSS descriptor is a
compact, batch-robust, scale-free summary of that multi-scale shape. No
test data and no test-time batch statistic are used: the rCLR is
strictly within-sample, and the scale ladder, kernels, standardisation
and head are all frozen here from training only.

## Usage

``` r
fit_inv_css(X_train, y_train, meta_train = NULL, hp = list())
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
  allowed-list): `n_scales` number of scales \\S\\ (default `6L`;
  integer \\\ge 2\\), `resample_len` fixed curve length \\L\\ (default
  `64L`; integer \\\ge 16\\; a very small \\L\\ relative to `sigma_max`
  thins the padding-free interior \\\[h+1, L-h-2\]\\ and can drive the
  largest scales into the neutral `(0,0,0,0)` descriptor – the default
  \\L=64\\ keeps every default scale well-populated), `sigma_min`
  smallest scale (default `1.0`; \\\>0\\), `sigma_max` largest scale
  (default `NULL` \\\to\\ `resample_len/8`; must exceed `sigma_min`),
  `shrink` ridge fraction of the mean covariance diagonal (default
  `0.1`; \\\ge 0\\), `min_features` per-specimen floor on the number of
  present rCLR coordinates (default `8L`; integer \\\ge 3\\), `eps`
  positive standard-deviation / ridge floor (default `1e-6`), and `seed`
  (default `1L`; accepted for interface uniformity only – the method
  consumes no randomness).

## Value

A plain list of class `inv_css_model` containing `feature_universe`, the
frozen scale ladder `scales` (per-scale `sigma`, half-width `h`,
normalised `kern`), `descriptor_dim` (\\= 4S\\), the frozen `head`
(`center`, `scale`, `active`, `w`, `b`, `ridge`), and the resolved `hp`.

## Details

The per-specimen score (see
[`score_inv_css`](https://kstawiski.github.io/OmicSelector/reference/score_inv_css.md))
is the frozen linear predictor on the standardised CSS descriptor,
\$\$S(x) = b + \sum\_{c} w_c \\ \frac{\varphi_c(x) - \mathrm{center}\_c}
{\mathrm{scale}\_c},\$\$ where \\\varphi(x)\\ is the specimen's CSS
descriptor. Every quantity except \\\varphi(x)\\ is frozen at fit time,
and \\\varphi(x)\\ depends only on \\x\\'s own rCLR profile, so the
score is computable from a single specimen and is exactly invariant to
per-sample positive scaling. Although it reads the same 1-D substrate as
the optimal-transport profile methods (the sorted rCLR curve), the
descriptor family here is a distinct multi-scale differential-geometric
(CSS) representation, not an optimal-transport embedding.

## References

Mokhtarian F, Mackworth AK. (1992) A theory of multiscale,
curvature-based shape representation for planar curves. *IEEE
Transactions on Pattern Analysis and Machine Intelligence* 14(8):
789-805.

Witkin AP. (1983) Scale-space filtering. *Proceedings of the 8th
International Joint Conference on Artificial Intelligence* 2: 1019-1022.

Martino C, Morton JT, Marotz CA, et al. (2019) A novel sparse
compositional technique reveals microbial perturbations. *mSystems*
4(1): e00016-19.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
p <- 60
features <- paste0("miR-", sprintf("%03d", seq_len(p)))
n <- 160
y <- rep(c(0, 1), each = n / 2)
X <- matrix(0, n, p, dimnames = list(NULL, features))
ctrl_mu <- seq(2, 5, length.out = p)                       # smooth ramp shape
case_mu <- ifelse(seq_len(p) <= p / 2, 2.75, 4.25)         # sharp step shape
for (i in seq_len(n)) {
  mu <- if (y[i] == 1) case_mu else ctrl_mu
  X[i, ] <- exp(stats::rnorm(p, mean = mu, sd = 0.2))
}
model <- fit_inv_css(X, y)
score <- score_inv_css(model, X)
} # }
```
