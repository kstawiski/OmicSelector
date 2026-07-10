# Spectral-normalized neural Gaussian process, distance-aware logit (SNGP)

Single-sample within-cohort discriminator (family N) implementing the
Spectral-normalized Neural Gaussian Process (SNGP) of Liu et al. (2020),
adapted to compositional specimens. A small **LayerNorm** MLP encoder
with **spectral-normalized** hidden Linear layers is trained on the
TRAINING set (via reticulate-python torch) to embed the per-sample-rCLR
specimen, then **FROZEN**; a FIXED random-Fourier-feature (RFF) map
approximates an RBF Gaussian process on the embedding, and a closed-form
LAPLACE posterior on the GP output weights gives the distance-aware
posterior-mean logit. A new specimen is scored by its
mean-field-adjusted GP logit. Larger score = more case-like.

**Distance awareness (Liu 2020).** SNGP combines two ingredients to make
a deterministic deep net distance-aware. (i) The hidden Linear layers
are **spectral-normalized** so the encoder is approximately bi-Lipschitz
and preserves input distances in the embedding space; combined with
**LayerNorm** (which normalises ACROSS FEATURES WITHIN A ROW, so it is
per-row and single-sample-exact – BatchNorm, which couples rows, is
forbidden here). (ii) A Gaussian-process output layer, approximated by
random Fourier features with a closed-form Laplace posterior, whose
predictive variance grows with distance from the training data, so
far-from-data specimens get a shrunken (less confident) logit. This is
exactly the SNGP recipe.

**Torch is confined to FIT.** The encoder is trained with torch (Adam,
CE loss) on the rCLR'd, universe-aligned training matrix with
`torch.nn.utils.parametrizations.spectral_norm` on the hidden Linear
layers. After training, the *effective* (post-spectral-norm,
coefficient- scaled) Linear weights \\W\_{\mathrm{eff}} =
c\\W\_{\mathrm{SN}}\\ – the exact weights torch used in its forward –
the biases, and the LayerNorm \\\gamma,\beta,\varepsilon\\ are EXPORTED
to R as plain numeric matrices/vectors; the classification head and the
python module are DISCARDED. The fixed RFF projection
\\W\_{\mathrm{rff}}, b\_{\mathrm{rff}}\\ (drawn ONCE under a dedicated
seeded torch generator) and the Laplace posterior (\\\beta\\,
\\\Sigma\\) are likewise frozen R-side. The fitted model therefore holds
NO external pointer and NO python dependency at score time: the encoder
forward (incl. LayerNorm), the RFF cosine map, and the mean-field GP
logit are all pure base-R. This makes the default `model_digest` a
stable deterministic snapshot and lets the package suite run the scorer
on any host without the venv.

**Spectral-norm bake-at-export fidelity.** PyTorch's `spectral_norm`
parametrisation materialises `layer$weight` to the spectrally-normalised
weight \\W\_{\mathrm{SN}}\\ (spectral norm \\1\\) at every forward; the
encoder's forward is \\x\\(c\\W\_{\mathrm{SN}})^\top + b\\ with the
bound coefficient \\c=\\`sn_coeff`. Since the training forward scales
the whole linear output by \\c\\ (\\c\\(H W\_{\mathrm{SN}}^\top + b) =
H\\(c W\_{\mathrm{SN}})^\top + c\\b\\), we export BOTH the effective
weight \\W\_{\mathrm{eff}} = c\\W\_{\mathrm{SN}}\\ and the effective
bias \\b\_{\mathrm{eff}} = c\\b\\, so the pure-R forward is a plain
matmul of the SAME weights AND bias torch used for ANY `sn_coeff` (not
only \\c=1\\): the R float64 forward matches the torch float32 forward
to \\\sim 10^{-7}\\ (the float32 round-off only).

**Fit/score consistency.** The Laplace posterior (\\\beta\\, \\\Sigma\\)
is computed from the PURE-R RFF features
\\\Phi(\\.\mathrm{unc\\sngp\\forward}(Z\_{\mathrm{train}})\\)\\ – the
IDENTICAL code path the scorer evaluates – NOT from a torch GP head. Fit
and score thus share one computation (the deepmaha / proto-net lesson:
fit and score must share one code path).

**The RFF Gaussian-process head.** With the frozen embedding \\z =
\mathrm{enc}(\mathrm{rclr}(x))\\ the random-Fourier feature map is
\$\$\Phi(z) = \sqrt{2/D}\\\cos(z\\W\_{\mathrm{rff}}^\top +
b\_{\mathrm{rff}}),\$\$ \\W\_{\mathrm{rff}}\sim\mathcal N(0, 1/\ell^2)\\
(length-scale \\\ell\\= `rff_ls`),
\\b\_{\mathrm{rff}}\sim\mathrm{Unif}(0, 2\pi)\\, \\D=\\`rff_dim`; this
approximates an RBF kernel GP. A ridge-regularised logistic regression
of the binary label on \\\Phi(Z\_{\mathrm{train}})\\ (a few Newton/IRLS
steps) gives the Laplace posterior MEAN \\\beta\\ and PRECISION
\\\Sigma^{-1} = \lambda I + \sum_i p_i(1-p_i)\\\phi_i\phi_i^\top\\
(\\\lambda=\\`ridge`), and \\\Sigma = (\Sigma^{-1})^{-1}\\.

**Distance-aware mean-field score.** The score of a new specimen is the
mean-field-approximated posterior logit (Liu 2020): \$\$s(x) =
\frac{\Phi(z)\cdot\beta} {\sqrt{1 +
(\pi/8)\\\Phi(z)^\top\Sigma\\\Phi(z)}}.\$\$ The denominator IS the
distance awareness: a specimen far from the training manifold has a
large RFF-GP predictive variance \\\Phi^\top\Sigma\Phi\\ and so a
shrunken (less confident) logit. We use the FULL mean-field denominator
(not the posterior mean alone): it is the defining SNGP predictive and
is the distance-aware mechanism; within-cohort ranking is essentially
preserved (Spearman \\\approx 1\\ vs the mean alone) while the score is
the faithful distance-aware logit. Larger = more case-like.

**Single-sample transform.** Each specimen is mapped to the
self-contained per-sample robust CLR over its OWN strictly-positive
support (geometric-mean centring on `v > 0`) in the FROZEN feature
universe; this is exactly invariant to per-specimen scaling and uses no
cross-row statistic. Universe features absent from a specimen carry the
neutral rCLR value `0`. The SAME rCLR transform is applied to the frozen
training rows at fit and to every query at score; combined with
LayerNorm (per-row), the fixed RFF map, and the frozen posterior, a
specimen's score depends only on that specimen and the frozen model –
single-row == batch EXACTLY, and it passes
[`singlesample_assert_row_equivariant`](https://kstawiski.github.io/OmicSelector/reference/singlesample_assert_row_equivariant.md).

**Degenerate-neutral.** A query returns the neutral score `0` if fewer
than `hp$min_features` universe features overlap the query columns (a
column-overlap floor, batch-independent), or it has empty positive
support over the frozen universe tested on the (pre-rCLR) aligned
ORIGINAL abundances (`!any(X_use[i, ] > 0)`). A FLAT all-equal-positive
composition maps to the rCLR origin but is a VALID specimen and is
scored normally (NOT floored) – the empty-support test is on the
ORIGINAL abundances, NOT on `all(z == 0)`.

## References

Liu JZ, Lin Z, Padhy S, Tran D, Bedrax-Weiss T, Lakshminarayanan B.
(2020) Simple and Principled Uncertainty Estimation with Deterministic
Deep Learning via Distance Awareness. *Advances in Neural Information
Processing Systems (NeurIPS)* 33. arXiv:2006.10108.
