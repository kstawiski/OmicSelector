# Fit unconstrained least-squares importance-fitting (uLSIF) density-ratio LRT

Learns a frozen two-class discriminator from training data by *directly*
estimating the class-conditional density ratio \\w(x) =
p\_{\mathrm{case}}(x) / p\_{\mathrm{control}}(x)\\ – without ever
estimating either density – using unconstrained Least-Squares Importance
Fitting (uLSIF; Kanamori, Hido & Sugiyama 2009). The estimated ratio is
itself the Neyman-Pearson-optimal class-conditional likelihood ratio, so
larger \\w(z)\\ means more case-like.

Each specimen is first mapped to a self-contained, per-sample robust
centred-log-ratio (rCLR) representation \\z\\, where \\z_j = \log v_j -
\mathrm{mean}\_{k:\\v_k\>0}\log v_k\\ on the nonzero parts of the
specimen and \\z_j = 0\\ on its zero parts. Because the centring uses
only that specimen's own nonzero values, \\z(c\\v) = z(v)\\ for any
\\c\>0\\: the representation – and therefore every downstream score – is
exactly invariant to per-sample scaling (library size / input amount).

The density ratio is modelled as a non-negative combination of Gaussian
basis functions \\\phi_l(z) = \exp(-\lVert z - c_l\rVert^2 /
(2\sigma^2))\\ centred at \\b\\ frozen kernel centres \\c_l\\. The
centres are a deterministic frozen subsample of the *case* rCLR anchors
(cap `hp$n_centers`); the bandwidth \\\sigma\\ is either supplied as
`hp$sigma` or learned by the median heuristic on the pairwise Euclidean
distances among the centres. Writing \\\Phi\_{\mathrm{control}}\\ and
\\\Phi\_{\mathrm{case}}\\ for the kernel matrices of the control / case
anchors against the centres, the closed-form uLSIF coefficients solve
the ridge-regularized normal equations \$\$(H + \lambda I)\\\alpha =
h,\qquad H = \tfrac{1}{n\_{\mathrm{control}}}
\Phi\_{\mathrm{control}}^\top \Phi\_{\mathrm{control}},\qquad h =
\mathrm{colMeans}(\Phi\_{\mathrm{case}}),\$\$ and the non-negative parts
\\\alpha \leftarrow \max(\alpha, 0)\\ are kept so that \\\hat w(z) =
\sum_l \alpha_l \phi_l(z) \ge 0\\ is a valid density ratio. There is no
iteration: \\H\\ is a Gram matrix (positive semidefinite), so \\H +
\lambda I\\ is positive definite for \\\lambda\>0\\ and the single
linear solve is exact (an increase-to-positive-definite ridge loop
guards only pathological numerical rank deficiency).

Fitting stores the raw case/control anchor rows and the raw centre rows
(so the kernel model can be re-derived consistently over any
present-feature subset at scoring), the frozen training feature
universe, the resolved hyperparameters, and – for inspection and tests –
the full-universe representation (`centers`, `sigma`, `alpha`). No test
data and no test-time batch statistic are used.

## Usage

``` r
fit_dre_ulsif(X_train, y_train, meta_train = NULL, hp = list())
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

  List of frozen hyperparameters: `sigma` Gaussian kernel bandwidth
  (default `NULL`, the median heuristic; otherwise a positive finite
  number); `lambda` positive uLSIF ridge (default `1e-3`); `n_centers`
  positive-integer cap on the number of case-anchor kernel centres
  (default `100L`); `max_anchors_per_class` positive-integer anchor cap
  per class (default `200L`); `min_features` positive-integer
  feature-overlap floor at scoring (default `3L`); `eps` positive ridge
  increment / bandwidth-and-score floor (default `1e-6`); and `seed`
  non-negative integer used only for deterministic anchor / centre
  subsampling while restoring the global RNG state (default `1L`).

## Value

A plain list of class `dre_ulsif_model` containing `case_anchors`,
`control_anchors`, `center_anchors`, `feature_universe`, `repr` (the
full-universe frozen centres/bandwidth/coefficients), and the resolved
`hp`.

## Details

The score (see
[`score_dre_ulsif`](https://kstawiski.github.io/OmicSelector/reference/score_dre_ulsif.md))
is the log density ratio \\S(z) = \log\max(\hat w(z), \varepsilon)\\
with \\\hat w(z) = \sum_l \alpha_l \exp(-\lVert z - c_l\rVert^2 /
(2\sigma^2))\\. The log scale is the canonical likelihood-ratio-test
statistic \\\log(p\_{\mathrm{case}}/p\_{\mathrm{control}})\\; it is a
strictly monotone transform of \\\hat w\\ (so it preserves the
density-ratio ranking and discrimination) but compresses the heavy right
tail of the ratio, which keeps the statistic numerically stable. The
\\\varepsilon\\ floor (which is only ever active when the non-negative
kernel model evaluates to (near) zero) keeps the logarithm finite.

Degenerate inputs stay finite without caller-side special-casing. A
feature that is constant across a class's anchors contributes a constant
rCLR column that simply shifts the kernel arguments; the Gram matrix
\\H\\ remains positive semidefinite and the ridge keeps \\H + \lambda
I\\ invertible. If the kernel centres collapse to a single point (or are
all identical), the median heuristic has no pairwise distance to use and
the bandwidth falls back to a within-centre root-mean-square length
scale floored at `hp$eps`, so \\\sigma \> 0\\ always.

## References

Kanamori T, Hido S, Sugiyama M. (2009) A least-squares approach to
direct importance estimation. *Journal of Machine Learning Research* 10:
1391-1445.

Sugiyama M, Suzuki T, Kanamori T. (2012) *Density Ratio Estimation in
Machine Learning*. Cambridge University Press.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
n <- 120; p <- 30; k <- 10
L <- matrix(stats::rnorm(n * p, 4, 0.5), nrow = n,
            dimnames = list(NULL, paste0("miR-", seq_len(p))))
y <- rep(c(0, 1), each = n / 2)
L[y == 1, seq_len(k)] <- L[y == 1, seq_len(k)] + 1.2  # case density shift
X <- exp(L)
model <- fit_dre_ulsif(X, y)
score <- score_dre_ulsif(model, X)
} # }
```
