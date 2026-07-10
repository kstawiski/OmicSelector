# Fit MFDFA spectrum (frac-mfdfa) within-sample discriminator

Learns a frozen, single-sample discriminator from the multifractal
detrended fluctuation analysis (MFDFA) spectrum of a specimen's
canonical-order robust-CLR (rCLR) abundance landscape. The construction
has six stages:

1.  **Frozen canonical feature order.** Features are sorted once by
    decreasing training column mean, with ties broken by feature index.
    The resulting `order_perm` is reused unchanged at scoring.

2.  **Per-specimen rCLR landscape.** For one specimen, model-universe
    features present in `X` are placed in the frozen order. Over that
    specimen's own positive support \\P=\\j:v_j\>0\\\\, the
    within-sample rCLR is \\g_j=\log v_j-\mathrm{mean}\_{k\in P}\log
    v_k\\. Structural zeros are excluded, but coordinates whose
    centred-log value is exactly zero are kept. This landscape is
    exactly invariant to multiplying the specimen by any positive
    scalar.

3.  **Fixed-length resample.** MFDFA needs one frozen absolute scale
    ladder. Because the positive support size can vary by specimen, the
    rCLR landscape is linearly interpolated to `hp$resample_len` equally
    spaced positions over the native index range. This imposes a fixed
    low-pass at the smallest scales, so the default ladder starts at
    `s_min >= 8` and avoids the smallest windows.

4.  **Base-R MFDFA.** The length-\\L\\ series is profiled as
    \\Y(i)=\sum\_{k\le i}(x_k-\bar x)\\. For each frozen scale,
    non-overlap windows are taken from both the start and the end of the
    profile. A polynomial trend of order `hp$m_poly` is removed by least
    squares, the window mean residual square is floored at `hp$eps`, and
    fluctuation functions \\F_q(s)\\ are formed for the frozen `hp$q`
    grid. The generalized Hurst exponent \\h(q)\\ is the slope of \\\log
    F_q(s)\\ against \\\log s\\.

5.  **Multifractal descriptors.** The fixed descriptor contains all
    \\h(q)\\ values, spectrum width `dAlpha`, \\h(2)\\, and a spectrum
    asymmetry term
    \\\alpha(q\_{max})+\alpha(q\_{min})-2\alpha((q\_{min}+q\_{max})/2)\\,
    with the midpoint alpha obtained by linear interpolation on the
    frozen q grid.

6.  **Frozen standardize + ridge-LDA head.** Descriptor center, scale,
    active flags, ridge, weight vector, and intercept are learned from
    training descriptors only. The head is oriented so larger scores are
    more case-like.

## Usage

``` r
fit_frac_mfdfa(X_train, y_train, meta_train = NULL, hp = list())
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

  Optional list of frozen MFDFA hyperparameters. Supported names are
  validated by the internal resolver: resample_len, n_scales, s_min,
  s_max, m_poly, q, shrink, min_features, eps, and seed.

## Value

A plain list of class `frac_mfdfa_model` containing `feature_universe`,
`order_perm`, frozen `scales`, `q`, `descriptor_dim`,
`descriptor_names`, frozen `head`, and resolved `hp`.

## Details

The per-specimen score is the frozen linear predictor \$\$S(x)=b+\sum_c
w_c\\\phi_c(x)-center_c\\/scale_c,\$\$ where \\\phi(x)\\ is the fixed
MFDFA descriptor. Every quantity except \\\phi(x)\\ is frozen at fit
time. At scoring, \\\phi(x)\\ is computed from only that row's present
positive values in the frozen feature order; no scored-batch statistic,
cross-row coupling, or test-batch normalization is used. The rCLR
landscape, fixed-length interpolation, MFDFA fluctuation ratios,
descriptors, and final score are invariant to per-sample positive
scaling.

## References

Kantelhardt JW, Zschiegner SA, Koscielny-Bunde E, et al. (2002)
Multifractal detrended fluctuation analysis of nonstationary time
series. *Physica A* 316(1-4): 87-114.

Peng CK, Buldyrev SV, Havlin S, et al. (1994) Mosaic organization of DNA
nucleotides. *Physical Review E* 49(2): 1685-1689.

Martino C, Morton JT, Marotz CA, et al. (2019) A novel sparse
compositional technique reveals microbial perturbations. *mSystems*
4(1): e00016-19.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
p <- 180
n <- 160
features <- paste0("miR-", sprintf("%03d", seq_len(p)))
y <- rep(c(0, 1), each = n / 2)
t <- seq(0, 1, length.out = p)
X <- matrix(0, n, p, dimnames = list(NULL, features))
for (i in seq_len(n)) {
  rough <- if (y[i] == 1) 0.8 * sin(48 * pi * t) else 0
  mu <- seq(5, 2, length.out = p) + rough
  X[i, ] <- exp(stats::rnorm(p, mu, 0.08))
}
model <- fit_frac_mfdfa(X, y)
score <- score_frac_mfdfa(model, X)
} # }
```
