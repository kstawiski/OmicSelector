# Fit topological persistence-image (tda-ph) within-sample discriminator

Learns a frozen, single-sample discriminator from the 0-dimensional
sublevel-set persistent homology of a specimen's within-sample rCLR
landscape. The construction has five stages, all estimated from training
data only (or, like the persistence-image Gaussian integral stencil,
fixed data-free constants):

1.  **Canonical feature order.** Features are arranged once, in
    *decreasing* training column mean of `X_train`; ties are broken by
    feature index, making `order_perm` deterministic and
    locale-independent. This fixed order defines the path graph on which
    every specimen's landscape is read.

2.  **Single-specimen rCLR landscape.** For one specimen, the
    model-universe features present in `X` are taken in the frozen order
    (`intersect(order_perm, colnames(X))`). Over that specimen's own
    positive support \\P=\\j: v_j\>0\\\\, the landscape is \\g_j=\log
    v_j-\mathrm{mean}\_{k\in P}\log v_k\\, keeping the frozen order and
    dropping only structural zeros. A coordinate with \\g_j=0\\ is kept;
    exclusion keys on \\v_j=0\\, never on \\g_j=0\\. This rCLR landscape
    is exactly invariant to per-sample positive scaling.

3.  **Path lower-star persistence.** The 0-dimensional finite
    sublevel-set classes of \\g\\ on the path graph are computed by a
    deterministic union-find merge tree. Vertices enter in increasing
    \\g\\; ties are broken by ascending position. When two existing
    components merge, the younger component (higher birth, ties by later
    birth position) dies at the saddle value. The essential
    global-minimum class is excluded by design.

4.  **Frozen persistence image.** Finite classes are mapped to
    birth-persistence coordinates and rasterised onto a frozen \\R\times
    R\\ grid by the exact Gaussian pixel integral using
    [`stats::pnorm`](https://rdrr.io/r/stats/Normal.html). The weight is
    \\w(p)=p\\. Vectorisation is row-major: pixel \\(a,c)\\ (birth pixel
    \\a\\, persistence pixel \\c\\) is stored at `(a - 1) * R + c`.

5.  **Frozen standardize + LDA head.** Training images are standardised
    by frozen per-pixel `center`/`scale`; exactly zero-variance pixels
    are inactive and contribute `0`. A ridge-regularised LDA head is fit
    on active standardised pixels, oriented so larger scores are more
    case-like.

The key topological design decision is that only *finite* 0-dimensional
classes are used. The essential global-minimum class is omitted because
it primarily re-encodes global rCLR spread; finite classes isolate
secondary valley structure in the canonical landscape. A strictly
monotone landscape has no finite classes and gives an all-zero
persistence image. That all-zero image is *not* short-circuited: it is
standardised and passed through the frozen head exactly as every other
descriptor in this family, so a monotone specimen receives the head's
learned, constant linear readout for the “no secondary valleys”
descriptor – a single deterministic baseline shared by every
empty-diagram specimen. This is the intended behaviour (the absence of
finite valleys is itself a discriminative cue, not an abstention); it is
*not* the bare intercept \\b\\, and it is distinct from the hard neutral
`0` that
[`score_tda_ph`](https://kstawiski.github.io/OmicSelector/reference/score_tda_ph.md)
returns only when a specimen falls below the `min_features` floor (a
genuine abstention on insufficient support). In the degenerate case
where *every* training specimen is monotone, all training images are
zero, every descriptor is deactivated, and the head collapses to the
constant-`0` scorer (`w = 0`, `b = 0`).

## Usage

``` r
fit_tda_ph(X_train, y_train, meta_train = NULL, hp = list())
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
  allowed-list): `grid_res` persistence-image side length \\R\\ (default
  `8L`; integer in `[3, 16]`), `sigma_frac` Gaussian bandwidth as a
  multiple of the smaller pixel width (default `1.0`; positive finite),
  `grid_pad` fractional padding of frozen grid ranges (default `0.05`;
  non-negative finite), `shrink` ridge fraction of the mean covariance
  diagonal (default `0.1`; non-negative finite), `min_features`
  per-specimen positive-support floor (default `8L`; integer \\\ge 3\\),
  `eps` positive standard-deviation / ridge floor (default `1e-6`), and
  `seed` (default `1L`; accepted for interface uniformity only – no
  randomness is used).

## Value

A plain list of class `tda_ph_model` containing `feature_universe`, the
frozen canonical `order_perm`, frozen persistence-image `grid`,
`descriptor_dim` (\\= R^2\\), frozen ridge-LDA `head` (`center`,
`scale`, `active`, `w`, `b`, `ridge`), and the resolved `hp`.

## Details

The per-specimen score (see
[`score_tda_ph`](https://kstawiski.github.io/OmicSelector/reference/score_tda_ph.md))
is the frozen linear predictor on the standardised persistence image,
\$\$S(x) = b + \sum_c w_c \frac{\phi_c(x)-\mathrm{center}\_c}
{\mathrm{scale}\_c},\$\$ where \\\phi(x)\\ is the specimen's
finite-class persistence image. Every quantity except \\\phi(x)\\ is
frozen at fit time, and \\\phi(x)\\ depends only on \\x\\'s own rCLR
landscape, its own 0-dimensional finite diagram, and the frozen image
grid. The score is therefore computable from one specimen alone, with no
test-batch normalisation, no cross-row coupling, and exact invariance to
per-sample positive scaling. A specimen whose finite (birth,
persistence) coordinates fall outside the frozen training grid has part
of its persistence-image mass truncated by the grid edges; the score is
still deterministic, finite, and rank/AUC-valid, but its raw magnitude
is not calibrated across that grid boundary.

## References

Edelsbrunner H, Letscher D, Zomorodian A. (2002) Topological persistence
and simplification. *Discrete & Computational Geometry* 28(4): 511-533.

Adams H, Emerson T, Kirby M, et al. (2017) Persistence images: a stable
vector representation of persistent homology. *Journal of Machine
Learning Research* 18(8): 1-35.

Martino C, Morton JT, Marotz CA, et al. (2019) A novel sparse
compositional technique reveals microbial perturbations. *mSystems*
4(1): e00016-19.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
p <- 72
n <- 160
features <- paste0("miR-", sprintf("%03d", seq_len(p)))
y <- rep(c(0, 1), each = n / 2)
X <- matrix(0, n, p, dimnames = list(NULL, features))
base_mu <- seq(10, 2, length.out = p)
block <- 24:48
pattern <- rep(c(0.45, -0.75, 0.35, -0.65), length.out = length(block))
for (i in seq_len(n)) {
  mu <- base_mu
  if (y[i] == 1) mu[block] <- mu[block] + pattern
  X[i, ] <- exp(stats::rnorm(p, mean = mu, sd = 0.02))
}
model <- fit_tda_ph(X, y)
score <- score_tda_ph(model, X)
} # }
```
