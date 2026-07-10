# Fit per-feature KDE naive-Bayes class-conditional LRT discriminator

Learns a frozen two-class discriminator from training data using a
*per-feature kernel density estimate* (KDE) of each class-conditional
marginal under a naive-Bayes (feature-independence) assumption. Each
specimen is first mapped to a self-contained, per-sample robust
centred-log-ratio (rCLR) representation \\z\\, where \\z_j = \log v_j -
\mathrm{mean}\_{k:\\v_k\>0}\log v_k\\ on the nonzero parts of the
specimen and \\z_j = 0\\ on its zero parts. Because the centring uses
only that specimen's own nonzero values, \\z(c\\v) = z(v)\\ for any
\\c\>0\\: the representation – and therefore every downstream score – is
exactly invariant to per-sample scaling (library size / input amount).

This is the *marginal* complement to the Gaussian-copula LRT. The copula
matches out the per-feature marginals via a probability-integral
transform and scores only the residual feature *dependence*; `lrt-nbkde`
does the opposite – it scores the per-feature class-conditional
*marginal* density ratio and ignores dependence (the naive-Bayes
assumption). For each feature \\j\\ and class \\c \in \\\mathrm{case},
\mathrm{control}\\\\ the frozen marginal is a Gaussian-kernel KDE of the
training rCLR values \\z^{\mathrm{train}}\_{c,j}\\: \$\$\log f\_{c,j}(t)
= \mathrm{logSumExp}\_{i}\\\left( -\tfrac{1}{2}\Big(\tfrac{t -
z^{\mathrm{train}}\_{c,j,i}}{h_j}\Big)^2\right) - \log n\_{c} - \log
h_j - \tfrac{1}{2}\log(2\pi),\$\$ evaluated directly at the specimen's
own coordinate \\t = z_j\\ (not on a fixed grid), so it is exact and
computable from a single specimen.

The per-feature bandwidth \\h_j\\ is frozen at fit. By default it is the
Silverman rule-of-thumb
[`stats::bw.nrd0`](https://rdrr.io/r/stats/bandwidth.html) of the
*pooled* rCLR values of feature \\j\\ across *both* classes (so the case
and control marginals of a feature are smoothed comparably), scaled by
`hp$bw_adjust` and floored at `hp$eps`: \\h_j =
\max(\mathrm{adjust}\cdot \mathrm{bw.nrd0}(\mathrm{pool}\_j),\\
\epsilon)\\. If `hp$bw` is supplied, that positive scalar replaces
`bw.nrd0` for every feature.

Fitting stores the raw case/control anchor rows (so the per-feature KDEs
can be re-derived consistently over any present-feature subset at
scoring), the frozen training feature universe, the resolved
hyperparameters, and – for inspection and tests – the full-universe
class representation (per-class rCLR anchor values and the per-feature
bandwidths). No test data and no test-time batch statistic are used.

## Usage

``` r
fit_lrt_nbkde(X_train, y_train, meta_train = NULL, hp = list())
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

  List of frozen hyperparameters: `bw` either `NULL` (default; the
  per-feature pooled
  [`stats::bw.nrd0`](https://rdrr.io/r/stats/bandwidth.html) bandwidth
  is used) or a single positive finite scalar bandwidth applied to every
  feature; `bw_adjust` positive finite multiplier on the bandwidth
  analogous to `density(adjust=)` (default `1.0`);
  `max_anchors_per_class` positive-integer anchor cap per class (default
  `200L`); `min_features` positive-integer feature-overlap floor at
  scoring (default `3L`); `eps` positive bandwidth/density floor
  (default `1e-6`); and `seed` non-negative integer used only for
  deterministic anchor subsampling while restoring the global RNG state
  (default `1L`).

## Value

A plain list of class `lrt_nbkde_model` containing `case_anchors`,
`control_anchors`, `feature_universe`, `repr` (the full-universe frozen
class representation: `z_case`, `z_control`, and the per-feature `bw`),
and the resolved `hp`.

## Details

The single-specimen score (see
[`score_lrt_nbkde`](https://kstawiski.github.io/OmicSelector/reference/score_lrt_nbkde.md))
is the naive-Bayes class-conditional marginal log-likelihood ratio
\$\$S(z) = \sum\_{j} \big\[\log f\_{\mathrm{case},j}(z_j) - \log
f\_{\mathrm{control},j}(z_j)\big\],\$\$ a sum of per-feature
class-conditional marginal log-density ratios; larger means more
case-like. The bandwidth and the \\-\tfrac12\log(2\pi)\\ constant cancel
within each feature's ratio, so only the kernel sums and the per-class
log-counts \\\log n_c\\ survive.

Degenerate features are handled without any caller-side special-casing.
A feature that is constant across a class's training anchors (e.g. zero
in all of them) still yields a finite KDE: the pooled `bw.nrd0` stays
positive (it falls back to a unit scale for a constant input) and is in
any case floored at `hp$eps`, and the Gaussian kernel sum is evaluated
through a stable log-sum-exp so each \\\log f\_{c,j}\\ is finite. The
fit and every score therefore remain finite.

## References

Silverman B.W. (1986) *Density Estimation for Statistics and Data
Analysis*. Chapman and Hall, London.

Scott D.W. (1992) *Multivariate Density Estimation: Theory, Practice,
and Visualization*. Wiley, New York.

Hand D.J., Yu K. (2001) Idiot's Bayes – not so stupid after all?
*International Statistical Review* 69(3): 385-398.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
n <- 120; p <- 30
L <- matrix(stats::rnorm(n * p, 4, 0.6), nrow = n,
            dimnames = list(NULL, paste0("miR-", seq_len(p))))
y <- rep(c(0, 1), each = n / 2)
L[y == 1, 1:8] <- L[y == 1, 1:8] + 1.2   # case-only marginal location shift
X <- exp(L)
model <- fit_lrt_nbkde(X, y)
score <- score_lrt_nbkde(model, X)
} # }
```
