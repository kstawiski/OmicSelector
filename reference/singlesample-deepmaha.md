# Deep-feature class-conditional Mahalanobis LRT (learned frozen embedding)

Single-sample within-cohort discriminator (family N) implementing the
class-conditional Mahalanobis confidence score of Lee et al. (2018),
adapted to compositional specimens. A small BatchNorm-free MLP embedding
of the per-sample-rCLR specimen is learned on the TRAINING set (via
reticulate-python torch), then **FROZEN**; per-class Gaussian means and
a single TIED covariance are estimated in the embedding space, and a new
specimen is scored by the class-conditional Mahalanobis log-likelihood
ratio in that frozen embedding. Larger score = more case-like.

**Torch is confined to FIT.** The MLP is trained with torch (Adam, CE
loss) on the rCLR'd, universe-aligned training matrix. After training,
the weight matrices and biases up to and INCLUDING the penultimate
embedding layer are EXPORTED to R as plain numeric matrices/vectors; the
final 2-logit classification head (not needed for scoring) and the
python module are DISCARDED. The fitted model therefore holds NO
external pointer and NO python dependency at score time: the embedding
forward, the Gaussian statistics, and the Mahalanobis LRT are all pure
base-R. This makes the default `model_digest` a stable deterministic
snapshot and lets the package suite run the scorer on any host without
the venv.

**Fit/score consistency despite float32.** Torch trains in float32; an R
float64 forward on the exported (float32-valued) weights differs from
the torch penultimate embedding by \\\sim 10^{-6}\\. To guarantee that
the frozen Gaussian statistics live in EXACTLY the same embedding the
scorer evaluates, the class means \\\mu\_{\mathrm{case}},
\mu\_{\mathrm{ctrl}}\\ and the tied covariance \\\Sigma\\ are computed
from the PURE-R forward of the training rows (`.lrt_deepmaha_forward`),
NOT from the torch embeddings. Fit and score thus use the identical
computation (the same lesson as nc-ecod-copod / nc-sinkhorn-single: fit
and score must share one code path).

**The embedding forward.** The exported net is \\\mathrm{input}(p) \to
\[\mathrm{Linear} \to \mathrm{ReLU}\]\_{1} \to \dots \to
\[\mathrm{Linear} \to \mathrm{ReLU}\]\_{L-1} \to \mathrm{Linear} \to
\mathrm{embedding}(d)\\. There are `length(hidden)` linear layers up to
the embedding (the last `hidden` width is the embedding dim \\d\\); the
activation (`"relu"` default, or `"tanh"`) is applied AFTER every linear
layer EXCEPT the final embedding layer (no activation follows the
embedding – it is the raw penultimate representation, matching the torch
net). `.lrt_deepmaha_forward(M, weights)` reimplements this exactly:
successive \\M \gets \phi(M W_l^\top + b_l)\\ for the hidden layers and
a final linear map to the \\d\\-dim embedding.

**Class-conditional Mahalanobis LRT.** With the tied pooled within-class
covariance \\\Sigma\\ (ridge-regularised) and class means
\\\mu\_{\mathrm{case}}, \mu\_{\mathrm{ctrl}}\\, the score of an embedded
specimen \\z\\ is \$\$s(z) = -\tfrac{1}{2}(z-\mu\_{\mathrm{case}})^\top
\Sigma^{-1} (z-\mu\_{\mathrm{case}}) +
\tfrac{1}{2}(z-\mu\_{\mathrm{ctrl}})^\top \Sigma^{-1}
(z-\mu\_{\mathrm{ctrl}}),\$\$ the class-conditional Gaussian
log-likelihood ratio (case minus control). With a single SHARED
covariance this is exactly LINEAR in \\z\\ (LDA in embedding space) –
the faithful Lee-2018 estimator, which also sidesteps \\p \gg n\\ since
\\d\\ is small. Larger = more case-like.

**Single-sample transform.** Each specimen is mapped to the
self-contained per-sample robust CLR over its OWN strictly-positive
support (geometric-mean centring on `v > 0`) in the FROZEN feature
universe; this is exactly invariant to per-specimen scaling and uses no
cross-row statistic. Universe features absent from a specimen carry the
neutral rCLR value `0`. The SAME rCLR transform is applied to the frozen
training rows at fit and to every query at score, so a specimen's score
depends only on that specimen and the frozen model – it passes
[`singlesample_assert_row_equivariant`](https://kstawiski.github.io/OmicSelector/reference/singlesample_assert_row_equivariant.md).

**Degenerate-neutral.** A query returns the neutral score `0` if fewer
than `hp$min_features` universe features overlap the query columns (a
column-overlap floor, batch-independent), or it has empty positive
support over the frozen universe tested on the (pre-rCLR) aligned
ORIGINAL abundances (`!any(X_use[i, ] > 0)`, where `X_use` is the query
aligned to the frozen universe with absent features set to 0). A FLAT
all-equal-positive composition maps to the rCLR origin but is a VALID
specimen and is scored normally (NOT floored) – the rCLR-origin trap
(the empty-support test is on the ORIGINAL abundances, NOT on
`all(z == 0)`).

## References

Lee K, Lee K, Lee H, Shin J. (2018) A Simple Unified Framework for
Detecting Out-of-Distribution Samples and Adversarial Attacks. *Advances
in Neural Information Processing Systems (NeurIPS)* 31.
arXiv:1807.03888.
