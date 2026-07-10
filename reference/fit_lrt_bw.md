# Fit a Bures-Wasserstein Gaussian-OT class-conditional LRT discriminator

Learns a frozen two-class discriminator from training data using a
Gaussian model per class whose covariance is regularized by
optimal-transport (Bures-Wasserstein, BW) geometry. Each specimen is
first mapped to a self-contained, per-sample robust centred-log-ratio
(rCLR) representation \\z\\, where \\z_j = \log v_j -
\mathrm{mean}\_{k:\\v_k\>0}\log v_k\\ on the nonzero parts of the
specimen and \\z_j = 0\\ on its zero parts. Because the centring uses
only that specimen's own nonzero values, \\z(c\\v) = z(v)\\ for any
\\c\>0\\: the representation – and therefore every downstream score – is
exactly (to numerical tolerance) invariant to per-sample scaling
(library size / input amount).

## Usage

``` r
fit_lrt_bw(X_train, y_train, meta_train = NULL, hp = list())
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

  List of frozen hyperparameters: `shrink` BW-geodesic shrink of each
  class covariance toward the shared anchor, in \\\[0,1)\\ (default
  `0.5`); `anchor` the shared anchor covariance, one of
  `"bw_barycenter"` or `"pooled"` (default `"bw_barycenter"`);
  `max_anchors_per_class` positive-integer anchor cap per class (default
  `200L`); `min_features` positive-integer feature-overlap floor at
  scoring (default `3L`); `eps` positive eigenvalue floor / ridge seed
  (default `1e-6`); and `seed` non-negative integer used only for
  deterministic anchor subsampling while restoring the global RNG state
  (default `1L`).

## Value

A plain list of class `lrt_bw_model` containing `case_anchors`,
`control_anchors`, `feature_universe`, `repr` (the full-universe frozen
class representation: per-class mean, regularized precision and log
determinant), and the resolved `hp`.

## Details

The single-specimen score (see
[`score_lrt_bw`](https://kstawiski.github.io/OmicSelector/reference/score_lrt_bw.md))
is the Gaussian log-likelihood ratio \$\$S(z) =
\ell\_{\mathrm{case}}(z) - \ell\_{\mathrm{control}}(z), \qquad \ell_c(z)
= -\tfrac{1}{2}(z-\mu_c)^\top \Sigma_c^{\mathrm{reg}\\-1}(z-\mu_c) -
\tfrac{1}{2}\log\det\Sigma_c^{\mathrm{reg}},\$\$ larger meaning more
case-like. The shared \\-\tfrac{p}{2}\log(2\pi)\\ term cancels in the
difference (both classes are scored over the same present feature set)
and is omitted.

Fitting stores the raw case/control anchor rows (so the class
representation can be re-derived consistently over any present-feature
subset at scoring), the frozen training feature universe, the resolved
hyperparameters, and – for inspection and tests – the full-universe
class representation. No test data and no test-time batch statistic are
used.

## Why a Gaussian likelihood and not a literal Wasserstein distance

The label *bw* here means **Bures-Wasserstein** – an optimal-transport
metric between Gaussians / SPD covariance matrices – and not a
bandwidth. A literal squared 2-Wasserstein distance from a *single*
specimen (a Dirac point, or a fixed-covariance Gaussian) to a class
Gaussian \\N(\mu_c,\Sigma_c)\\ is \\\lVert z - \mu_c\rVert^2 +
(\text{class-constant Bures term})\\: the class covariance enters *only*
through a per-class additive constant, so the W2-distance log-ratio
collapses to a Euclidean nearest-centroid (linear) classifier that does
*not* use \\\Sigma_c\\ discriminatively. To build a single-sample LRT
that genuinely uses the frozen class covariances we score with the
Gaussian log-likelihood ratio (Neyman-Pearson optimal for Gaussian class
models); the **Bures-Wasserstein element is carried by how the
covariances are estimated and regularized** – by OT-geometric
(BW-geodesic) shrinkage of each class covariance toward a shared
barycenter anchor. That is what distinguishes this method from a plain
Euclidean-shrinkage Mahalanobis LRT.

## Fit (training only)

On the per-sample rCLR \\z\\ of the frozen class anchors:

1.  \\\mu_c\\ is the class mean (the W2 / Frechet barycenter of the
    class points equals the ordinary Euclidean mean in rCLR space) and
    \\\Sigma_c\\ is the class sample covariance
    ([`stats::cov`](https://rdrr.io/r/stats/cor.html)), symmetrized.

2.  A shared anchor covariance \\\Sigma\_{\mathrm{pool}}\\ is either the
    equal-weight **Bures-Wasserstein barycenter** of
    \\\\\Sigma\_{\mathrm{case}},\Sigma\_{\mathrm{control}}\\\\ (the
    midpoint of the BW geodesic; `hp$anchor = "bw_barycenter"`) or the
    pooled within-class covariance (`hp$anchor = "pooled"`).

3.  Each class covariance is BW-geodesically shrunk toward the anchor by
    \\t = \\ `hp$shrink` \\\in \[0,1)\\: \\\Sigma_c^{\mathrm{reg}} =
    \mathrm{BWgeo}(\Sigma_c, \Sigma\_{\mathrm{pool}}, t)\\ (\\t=0\\
    keeps \\\Sigma_c\\; \\t=1\\ gives the anchor), then a diagonal ridge
    is added until `chol` succeeds. The Cholesky factor yields
    \\\Sigma_c^{\mathrm{reg}\\-1}\\ via `chol2inv` and
    \\\log\det\Sigma_c^{\mathrm{reg}} = 2\sum\log \mathrm{diag}\\.

## Bures-Wasserstein geodesic

For SPD matrices \\A,B\\ and \\t\in\[0,1\]\\ the McCann / OT
interpolation is \$\$T = A^{-1/2}\\(A^{1/2} B A^{1/2})^{1/2}\\A^{-1/2},
\quad M_t = (1-t) I + t\\T, \quad \mathrm{BWgeo}(A,B,t) =
M_t\\A\\M_t,\$\$ symmetrized. Symmetric matrix square roots use the
eigendecomposition \\A^{1/2} =
V\\\mathrm{diag}(\sqrt{\max(\lambda,0)})\\V^\top\\ and \\A^{-1/2} =
V\\\mathrm{diag}(1/\sqrt{\max(\lambda,\varepsilon)})\\V^\top\\ (tiny
eigenvalues floored at `hp$eps` for the inverse square root). Both \\A\\
and \\B\\ are floored to positive definite before forming the geodesic.
The equal-weight BW barycenter is \\\mathrm{BWgeo}(A,B,0.5)\\.

## High-dimension (p \> n) behaviour

When the number of features \\p\\ exceeds a class's anchor count
\\n_c\\, \\\Sigma_c\\ is rank-deficient. The eigenvalue floor
(null-space directions receive variance `hp$eps`), the BW-geodesic
shrinkage toward the better-conditioned anchor (with a substantial
default `hp$shrink = 0.5`), and the final ridge-to-positive-definite
loop together guarantee a positive definite \\\Sigma_c^{\mathrm{reg}}\\,
so the fit and every score remain finite.

## References

Bhatia R., Jain T., Lim Y. (2019) On the Bures-Wasserstein distance
between positive definite matrices. *Expositiones Mathematicae* 37(2):
165-191.

Agueh M., Carlier G. (2011) Barycenters in the Wasserstein space. *SIAM
Journal on Mathematical Analysis* 43(2): 904-924.

Takatsu A. (2011) Wasserstein geometry of Gaussian measures. *Osaka
Journal of Mathematics* 48(4): 1005-1026.

Neyman J., Pearson E.S. (1933) On the problem of the most efficient
tests of statistical hypotheses. *Philosophical Transactions of the
Royal Society A* 231: 289-337.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
n <- 160; p <- 30
L <- matrix(stats::rnorm(n * p, 4, 0.6), nrow = n,
            dimnames = list(NULL, paste0("miR-", seq_len(p))))
y <- rep(c(0, 1), each = n / 2)
shift <- rep(c(1.2, -1.2), length.out = 12)
L[y == 1, 1:12] <- L[y == 1, 1:12] + matrix(shift, sum(y == 1), 12, byrow = TRUE)
X <- exp(L)
model <- fit_lrt_bw(X, y)
score <- score_lrt_bw(model, X)
} # }
```
