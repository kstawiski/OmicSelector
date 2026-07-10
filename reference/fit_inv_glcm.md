# Fit GLCM Haralick texture (inv-glcm) within-sample discriminator

Learns a frozen, single-sample discriminator from the 1-D Haralick
texture of a specimen's canonical-order robust-CLR (rCLR) abundance
landscape. The construction has five stages:

1.  **Frozen canonical feature order.** Features are sorted once by
    decreasing training column mean (`colMeans(X_train)`), ties broken
    by feature index (locale-independent, deterministic). The resulting
    permutation `order_perm` is reused unchanged at scoring.

2.  **Per-specimen rCLR landscape.** Over each specimen's own positive
    support \\\\j : v_j \> 0\\\\, the within-sample robust-CLR is \\g_j
    = \log v_j - \mathrm{mean}\_{k \in P}\log v_k\\. Structural zeros
    (\\v_j = 0\\) are excluded; present coordinates with \\g_j = 0\\ are
    kept (keyed on `v > 0`, NOT `g != 0`). Values are placed in the
    frozen canonical order, yielding a sequence \\g_1, \dots, g_m\\. The
    rCLR landscape \\g\\ is exactly scale-invariant because \\\log(c
    v_j) - \mathrm{mean}\log(c v) = \log v_j - \mathrm{mean}\log v\\ for
    any \\c \> 0\\; the frozen-edge quantised levels, hence the score,
    are therefore unchanged (bit-identical in practice – the only
    exception is the measure-zero set where an rCLR value lies within
    floating-point rounding of a frozen edge and a positive rescaling
    could flip its bin).

3.  **Frozen gray-level quantisation.** `K = hp$n_levels` gray levels
    are defined by \\K-1\\ interior quantile edges learned from the
    pooled training landscape values. Each specimen's landscape is
    quantised via left-closed binning: `findInterval(g, edges) + 1`
    (intervals \\\[e_i, e\_{i+1})\\, so a value exactly equal to an edge
    falls into the higher-numbered level), yielding levels in \\\\1,
    \dots, K\\\\. Because both \\g\\ and the frozen `edges` live on the
    scale-invariant rCLR scale, the quantised level sequence is
    unchanged by positive per-sample scaling.

4.  **Per-specimen 1-D GLCM + Haralick descriptors.** From the quantised
    sequence \\q_1, \dots, q_m\\ a distance-1 co-occurrence matrix
    \\C\[i,j\]\\ is built (\\C\[q_t, q\_{t+1}\]\\ incremented for \\t =
    1, \dots, m-1\\), symmetrised (\\C \leftarrow C + C^\top\\), and
    normalised to a probability matrix \\P = C / \mathrm{sum}(C)\\. Five
    Haralick features are extracted: contrast, homogeneity (IDM), energy
    (ASM), entropy, and correlation (with degeneracy guard). This helper
    is defined once and used identically at fit and at score.

5.  **Frozen standardize + ridge-LDA head.** Per-descriptor frozen
    `center`/`scale` (sd floored at `eps`; zero-variance descriptors
    deactivated), pooled within-class covariance ridged to PD, \\w =
    S_w^{-1}(\mu\_\text{case} - \mu\_\text{ctrl})\\, \\b = -\frac12
    (\mu\_\text{case} + \mu\_\text{ctrl})^\top w\\. Larger score = more
    case-like.

## Usage

``` r
fit_inv_glcm(X_train, y_train, meta_train = NULL, hp = list())
```

## Arguments

- X_train:

  Numeric matrix (samples \\\times\\ features) of non-negative abundance
  values. Columns must be uniquely named feature ids.

- y_train:

  Integer/numeric 0/1 labels aligned to `X_train`; 1 is case/disease, 0
  is control.

- meta_train:

  Optional per-sample metadata. Accepted for interface uniformity;
  ignored.

- hp:

  List of frozen hyperparameters (strict allowed-list, exact `[[ ]]`
  reads): `n_levels` K (default `8L`; integer in `[2, 32]`) – GLCM gray
  levels; `shrink` (default `0.1`; non-negative finite) – ridge fraction
  of mean covariance diagonal; `min_features` (default `8L`; integer
  \\\ge 3\\) – per-specimen positive-support floor; `eps` (default
  `1e-6`; positive finite) – sd / ridge / correlation-degeneracy floor;
  `seed` (default `1L`; non-negative integer) – accepted for interface
  uniformity, no randomness is used.

## Value

A plain list of class `inv_glcm_model` containing `feature_universe`,
`order_perm` (feature names in frozen canonical order), `edges` (frozen
\\K-1\\ quantile breakpoints), `K` (number of gray levels),
`descriptor_dim` (5), the frozen `head` (`center`, `scale`, `active`,
`w`, `b`, `ridge`), and the resolved `hp`.

## Details

The per-specimen score is the frozen linear predictor on the
standardised Haralick descriptor: \$\$S(x) = b + \sum_c w_c
\frac{\phi_c(x) - \mathrm{center}\_c} {\mathrm{scale}\_c},\$\$ where
\\\phi(x)\\ is the 5-dim descriptor. Every quantity except \\\phi(x)\\
is frozen at fit time; \\\phi(x)\\ depends only on \\x\\'s own rCLR
landscape, so the score is computable from a single specimen and is
unchanged by per-sample positive scaling (bit-identical in practice; see
[`score_inv_glcm`](https://kstawiski.github.io/OmicSelector/reference/score_inv_glcm.md)
for the measure-zero quantiser-edge qualification).

## References

Haralick RM, Shanmugam K, Dinstein I. (1973) Textural features for image
classification. *IEEE Transactions on Systems, Man, and Cybernetics*
SMC-3(6): 610-621.

Martino C, Morton JT, Marotz CA, et al. (2019) A novel sparse
compositional technique reveals microbial perturbations. *mSystems*
4(1): e00016-19.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
p <- 50
X <- matrix(stats::rgamma(120 * p, shape = 20, rate = 0.2), nrow = 120,
            dimnames = list(NULL, paste0("miR-", sprintf("%03d", seq_len(p)))))
y <- rep(c(0, 1), each = 60)
X[y == 1, 11:20] <- X[y == 1, 11:20] + 150
model <- fit_inv_glcm(X, y)
score <- score_inv_glcm(model, X)
} # }
```
