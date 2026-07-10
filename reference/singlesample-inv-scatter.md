# 1D wavelet-scattering single-sample discriminator (python-at-score, frozen head)

Single-sample within-cohort discriminator (family M, image/invariance)
that maps each specimen's frozen-order per-sample robust-CLR (rCLR)
profile to a **1D wavelet scattering transform** (Mallat 2012; Anden &
Mallat 2014) and scores the scattering coefficients through a FROZEN
regularized-logistic head. A specimen of \\p\\ universe features is laid
out as a 1D signal along the frozen feature axis, zero-padded to
\\\mathrm{shape} = \mathrm{next\\pow2}(p)\\, scattered by a fixed
cascade of wavelet convolutions + modulus + low-pass averaging
(`ScatteringTorch1D`, kymatio), flattened to a fixed-length coefficient
vector, standardized by FROZEN training statistics, and fed to a FROZEN
ridge-logistic linear head. Larger = more case-like.

**Python-at-score contract (no live pointer).** Like the proven
python-at-score siblings
[`fit_img_gasfcnn`](https://kstawiski.github.io/OmicSelector/reference/fit_img_gasfcnn.md)
/
[`fit_tabicl`](https://kstawiski.github.io/OmicSelector/reference/fit_tabicl.md)
(and unlike the export-to-pure-R trained-torch siblings `proto-net` /
`lrt-deepmaha`), the scattering transform is computed in python at fit
AND score, but the fitted model holds ONLY plain R-side state: the
scattering CONFIG (`J`, `Q`, `shape`, `max_order`, the frozen feature
universe), the FROZEN per-coordinate coefficient standardization
(`coef_mean`, `coef_sd`), the FROZEN linear head (`weights`,
`intercept`), and the seed/hp. There is NO live python pointer on the
model. At SCORE the identical `ScatteringTorch1D` is REBUILT in python
from the config each call, the coefficients are computed FORCED
ROW-BY-ROW in float64, standardized with the frozen stats, and the
frozen R-side linear head is applied in pure base R. Because every model
field is a plain R object, the default `model_digest` is a stable
deterministic snapshot suitable for
[`singlesample_assert_row_equivariant`](https://kstawiski.github.io/OmicSelector/reference/singlesample_assert_row_equivariant.md).

**Why this is inherently single-sample (no BatchNorm).** The scattering
transform \\S(x)\\ is a FIXED cascade of wavelet convolutions, modulus,
and low-pass averaging with NO learned parameters and NO cross-sample
statistics – each specimen's coefficients depend ONLY on that specimen's
own signal. There is no BatchNorm and no batch coupling. The only
numerical wrinkle is float32 convolution accumulation across a batch
dimension; this scorer computes the coefficients in **float64** ONE ROW
AT A TIME, so the \\n=1\\ path is used uniformly and a row's
coefficients (hence its score) are bit-identical whether scored alone or
in any batch – single-row score == batch score EXACTLY 0. The
coefficient standardization uses FROZEN training mean/sd (never the
query's own), and the head is a frozen linear map, so the whole pipeline
is single-sample-safe.

**Signal construction (deterministic, frozen).** For a specimen's
universe-aligned per-sample rCLR vector \\v\\ (length \\p = \\ size of
the FROZEN feature universe, in the FROZEN feature ORDER), the 1D signal
is \\v\\ laid out along the frozen feature axis, zero-padded to
\\\mathrm{shape} = \mathrm{next\\pow2}(p)\\ (FROZEN at fit from the
training universe size). The scattering transform is ORDER-dependent
along the feature axis: this is an intended, frozen design choice (the
feature order is fixed at fit; there is NO claim of feature-permutation
invariance). The row-equivariance harness permutes SAMPLES (rows), never
features (columns) – exactly as for img-gasfcnn's GASF.

**Single-sample transform.** Each specimen is mapped to the
self-contained per-sample robust CLR over its OWN strictly-positive
support (geometric-mean centring on `v > 0`) in the FROZEN feature
universe; this is exactly invariant to per-specimen positive scaling and
uses no cross-row statistic. Universe features absent from a specimen
carry the neutral rCLR value `0`. The SAME rCLR transform feeds the
scattering at fit and at score.

**Degenerate-neutral.** A query returns the neutral score `0` if fewer
than `hp$min_features` universe features overlap the query columns (a
column-overlap floor, batch-independent), or it has empty positive
support over the frozen universe on the (pre-rCLR) aligned ORIGINAL
abundances (`!any(X_use[i, ] > 0)`). A FLAT all-equal-positive
composition maps to the rCLR origin (an all-zero signal) but is a VALID
specimen: it yields valid (zero) scattering coefficients and is scored
normally (NOT floored) – the empty-support test is on the ORIGINAL
abundances, NOT on `all(z == 0)`.

**kymatio import workaround.** kymatio 0.3.0's package init
(`import kymatio.torch`) eagerly imports the 3D harmonic-scattering
module, which calls `scipy.special.sph_harm` (REMOVED in scipy \>=
1.15), so the top-level import RAISES `ImportError`. The 1D frontend is
therefore imported DIRECTLY from
`kymatio.scattering1d.frontend.torch_frontend` (verified working:
float64, deterministic). The reticulate skip-guard probes `torch` AND
that direct `ScatteringTorch1D` import (NOT `import kymatio.torch`).

## References

Mallat S. (2012) Group Invariant Scattering. *Communications on Pure and
Applied Mathematics* 65(10):1331-1398.

Anden J, Mallat S. (2014) Deep Scattering Spectrum. *IEEE Transactions
on Signal Processing* 62(16):4114-4128.

Andreux M, Angles T, Exarchakis G, et al. (2020) Kymatio: Scattering
Transforms in Python. *Journal of Machine Learning Research* 21(60):1-6.
