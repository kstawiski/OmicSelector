# Fit a linearized-optimal-transport (CDT tangent) LRT discriminator

Learns a frozen two-class discriminator from training data by embedding
each specimen into **linearized-optimal-transport (LOT) tangent space**
– in one dimension the exact **Cumulative Distribution Transform (CDT)**
of Park-Kolouri-Rohde – and then separating the two classes with a
linear (LDA-style, pooled-covariance) Gaussian discriminant in that
tangent space. Each specimen is first mapped to a self-contained,
per-sample robust centred-log-ratio (rCLR) representation \\z\\, where
\\z_j = \log v_j - \mathrm{mean}\_{k:\\v_k\>0}\log v_k\\ on its nonzero
parts and \\z_j = 0\\ on its zero parts. Because the centring uses only
that specimen's own nonzero values, \\z(c\\v) = z(v)\\ for any \\c\>0\\:
the representation – and therefore every downstream score – is exactly
(to numerical tolerance) invariant to per-sample scaling (library size /
input amount).

## Usage

``` r
fit_ot_lot(X_train, y_train, meta_train = NULL, hp = list())
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

  List of frozen hyperparameters: `n_quantiles` CDT grid size \\K\\
  (positive integer \\\ge 4\\, default `16L`); `reference` the tangent
  base point, one of `"barycenter"` or `"pooled"` (default
  `"barycenter"`); `shrink` fixed diagonal-load fraction added as
  `shrink * mean(diag(S))` before the ridge loop (non-negative finite,
  default `0.1`); `min_features` minimum present rCLR coordinates needed
  to form a quantile function (positive integer \\\ge 2\\, default
  `5L`); `eps` positive ridge increment / seed (default `1e-6`); and
  `seed` non-negative integer kept for API consistency only – no RNG is
  consumed, and a fit leaves `.Random.seed` untouched (default `1L`).

## Value

A plain list of class `ot_lot_model` containing `t_k` (the frozen
quantile grid), `Q_ref` (the frozen reference quantile function), `head`
(the frozen linear discriminant: `w`, `b`, per-class tangent means
`mu_case`/`mu_control`, `Sinv`, the resolved ridge, and the per-class
usable counts), `feature_universe`, `reference`, the per-class usable
counts, and the resolved `hp`.

## Details

The single-specimen score (see
[`score_ot_lot`](https://kstawiski.github.io/OmicSelector/reference/score_ot_lot.md))
is the frozen linear functional \\S(Q) = w^\top(Q - Q\_{\mathrm{ref}}) +
b\\ of the specimen's own CDT embedding, larger meaning more case-like.
The reference \\Q\_{\mathrm{ref}}\\, the grid \\t_k\\ and the head
\\(w,b)\\ are all frozen from training; the specimen's quantile function
is computed from its own present coordinates, so no scored-batch
statistic and no cross-row coupling enter. Because the whole model lives
in the fixed-dimension CDT tangent space \\\mathbb{R}^K\\ (no
feature-indexed component beyond which features are admitted), partial
feature overlap introduces no systematic batch-statistic drift – there
is no cross-sample renormalisation, and the same specimen scored alone
always yields the same score (the inclusion-gate property). Fewer
present features do, however, lower the resolution of (and add sampling
variance to) the empirical quantile function, which is still evaluated
at the frozen \\t_k\\; this added variance is intrinsic to the
value-distribution representation, not a leak of the present-feature
count.

## The single-specimen measure (why a 1-D CDT and not a literal W2)

A literal squared 2-Wasserstein distance from a *single* specimen (a
Dirac point measure) to a class distribution is \\\lVert\cdot\rVert^2 +
\text{const}\\ – a covariance-blind nearest-centroid that discards the
class shape, the same collapse that forces a Gaussian LRT in
[`fit_lrt_bw`](https://kstawiski.github.io/OmicSelector/reference/fit_lrt_bw.md).
We avoid it by treating each specimen as a **one-dimensional empirical
measure over its own present rCLR coordinate values**. Take the present
rCLR coordinates \\z\_{\mathrm{pos}} = \\z_j: v_j \> 0\\\\ (rCLR sets
absent parts to exactly \\0\\; including the structural-zero point mass
would leak detection rate / sequencing depth – a batch artefact – so
only structural zeros, keyed on the raw abundance \\v_j = 0\\, are
excluded; a present coordinate that happens to equal the geometric mean
(\\z_j = 0\\) is kept). \\z\_{\mathrm{pos}}\\ is mean-zero by rCLR
construction, so its *location* is fixed at \\0\\ and the discriminative
content is the **distributional shape** (spread, skew, tails) of the
centred-log-ratio profile. The specimen's measure is \\\nu =
\tfrac{1}{m}\sum\_{j}\delta\_{z\_{\mathrm{pos},j}}\\ (\\m =
\|z\_{\mathrm{pos}}\|\\); its quantile function \\Q\_\nu\\ evaluated at
fixed interior levels by the empirical inverse-CDF (`type = 1`) is the
exact 1-D optimal-transport (monotone rearrangement) representation.
This representation is **feature-exchangeable by design** – it scores
the shape of the rCLR value-distribution – a mechanism distinct from and
complementary to the feature-identity-preserving sliced / Gaussian LRTs.

## Fit (training only; no RNG consumed)

1.  A fixed interior quantile grid \\t_k = (k-\tfrac12)/K\\, \\k =
    1,\dots,K\\ (\\K = \\ `hp$n_quantiles`).

2.  For each training specimen the per-sample rCLR is formed over the
    feature universe, its present coordinates \\z\_{\mathrm{pos}}\\ are
    taken, and – if at least `hp$min_features` of them exist – its CDT
    quantile function \\Q = \mathrm{quantile}(z\_{\mathrm{pos}}, t_k,
    \mathrm{type}=1)\\ (the empirical inverse-CDF, i.e. the exact 1-D OT
    map) \\\in \mathbb{R}^K\\ is recorded (a specimen with fewer usable
    coordinates is skipped; at least two usable specimens per class are
    required).

3.  The **frozen reference quantile function** \\Q\_{\mathrm{ref}}\\
    (the LOT tangent base point) is the 1-D Wasserstein barycenter
    \\\mathrm{colMeans}(Q_i)\\ over all usable specimens
    (`hp$reference = "barycenter"`; the 1-D barycenter quantile function
    is the mean of the quantile functions) or the pooled quantile
    function \\\mathrm{quantile}(\bigcup_i z\_{\mathrm{pos},i}, t_k)\\
    (`hp$reference = "pooled"`).

4.  The **CDT tangent embedding** is \\v_i = Q_i - Q\_{\mathrm{ref}}\\
    (the Wasserstein log-map at the reference).

5.  A **linear discriminant head** (LDA-style, which keeps the method
    faithful to *linear* OT): per-class mean \\\mu_c\\, pooled
    within-class covariance \\S = ((n_1-1)\\\mathrm{Cov}\_1 +
    (n_0-1)\\\mathrm{Cov}\_0)/(n_1+n_0-2)\\. Because the quantile vector
    is monotone, adjacent components are highly correlated and \\S\\ is
    near-singular, so regularization is required: a fixed diagonal load
    `hp$shrink`\\\times\mathrm{mean}(\mathrm{diag} S)\\ is added first,
    then a ridge is raised from `hp$eps` (doubling) until `chol`
    succeeds. With \\S\_{\mathrm{inv}} =
    \mathrm{chol2inv}(\mathrm{chol}(S+\lambda I))\\ the head is \\w =
    S\_{\mathrm{inv}}(\mu_1-\mu_0)\\, \\b = -\tfrac12(\mu_1+\mu_0)^\top
    w\\.

## Reference-origin invariance

The decision value \\w^\top(Q-Q\_{\mathrm{ref}})+b\\ is algebraically
independent of the choice of \\Q\_{\mathrm{ref}}\\: \\w\\ depends only
on the class-mean difference (a constant shift cancels) and the
\\Q\_{\mathrm{ref}}\\ terms in \\w^\top(Q-Q\_{\mathrm{ref}})\\ and \\b\\
cancel, leaving the textbook LDA decision \\w^\top Q - \tfrac12
w^\top(\bar Q_1+\bar Q_0)\\. The reference fixes the tangent origin (and
is reported for inspection) but does not change any score.

## References

Park S.R., Kolouri S., Kundu S., Rohde G.K. (2018) The cumulative
distribution transform and linear pattern classification. *Applied and
Computational Harmonic Analysis* 45(3): 616-641.

Kolouri S., Park S.R., Thorpe M., Slepcev D., Rohde G.K. (2017) Optimal
mass transport: signal processing and machine-learning applications.
*IEEE Signal Processing Magazine* 34(4): 43-59.

Wang W., Slepcev D., Basu S., Ozolek J.A., Rohde G.K. (2013) A linear
optimal transportation framework for quantifying and visualizing
variation in sets of images. *International Journal of Computer Vision*
101(2): 254-269.

Fisher R.A. (1936) The use of multiple measurements in taxonomic
problems. *Annals of Eugenics* 7(2): 179-188.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
n <- 200; p <- 40
base <- 6
y <- rep(c(0, 1), each = n / 2)
E <- matrix(0, n, p, dimnames = list(NULL, paste0("miR-", seq_len(p))))
E[y == 0, ] <- stats::rnorm(sum(y == 0) * p, 0, 0.30)   # tight controls
E[y == 1, ] <- stats::rnorm(sum(y == 1) * p, 0, 0.60)   # dispersed cases
X <- exp(base + E)
model <- fit_ot_lot(X, y)
score <- score_ot_lot(model, X)
} # }
```
