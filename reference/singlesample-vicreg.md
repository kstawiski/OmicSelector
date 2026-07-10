# VICReg self-supervised embedding + frozen linear-probe discriminator

Single-sample within-cohort discriminator (family N) implementing a
self-supervised VICReg representation in the sense of Bardes, Ponce &
LeCun (2022), adapted to compositional specimens. A small BatchNorm-free
MLP encoder of the per-sample-rCLR specimen is learned on the TRAINING
set by SELF-SUPERVISED **VICReg** regularization
(variance-invariance-covariance loss on two augmented views, **no
labels**) via reticulate-python torch, then **FROZEN**; a frozen
linear-probe head (the SSL linear-probe evaluation protocol) is fit on
the labeled training embeddings, and a new specimen is scored by its
frozen-head logit on its pure-R encoder embedding. Larger score = more
case-like.

**VICReg self-supervised training (Bardes 2022).** Unlike the
cross-entropy-trained `lrt-deepmaha` sibling or the prototypical-loss
`proto-net` sibling, the encoder here is trained with NO labels by the
VICReg loss. Each step, for a batch of (rCLR'd, universe-aligned)
training rows, TWO augmented views are formed (random feature masking +
additive Gaussian jitter in rCLR space, independent draws per view),
embedded through the encoder and a small expander/projector to \\Z,
Z'\\, and the VICReg loss (Bardes 2022 eq. 6) \$\$\mathcal{L} =
\lambda\\s(Z,Z') + \mu\\\[v(Z)+v(Z')\] + \nu\\\[c(Z)+c(Z')\]\$\$ is
minimised, where the invariance term \\s\\ is the MSE between \\Z\\ and
\\Z'\\; the variance term \\v(Z) = \frac1d \sum_j \max(0, \gamma -
\sqrt{\mathrm{Var}(Z_j) + \epsilon})\\ is the hinge (at \\\gamma = 1\\)
on the per-dimension embedding std; and the covariance term \\c(Z) =
\frac1d \sum\_{i\neq j} \[\mathrm{Cov}(Z)\_{ij}\]^2\\ penalises the
squared off-diagonal of the embedding covariance. Defaults \\\lambda =
25, \mu = 25, \nu = 1\\ (the paper's standard). Adam optimises encoder +
expander over `epochs` steps. The DOWNSTREAM embedding is the ENCODER
output (pre-expander) — standard VICReg practice; the expander is
discarded and only the encoder is exported.

**Torch is confined to FIT.** The encoder is trained with torch on the
rCLR'd, universe-aligned training matrix. After SSL training, the
encoder weight matrices and biases up to and INCLUDING the embedding
layer are EXPORTED to R as plain numeric matrices/vectors; the expander,
the python module, and the head's torch object are all DISCARDED (only
the head's exported \\w, b\\ survive). The fitted model therefore holds
NO external pointer and NO python dependency at score time: the encoder
forward and the linear-head logit are pure base-R. This makes the
default `model_digest` a stable deterministic snapshot and lets the
package suite run the scorer on any host without the venv.

**Fit/score consistency despite float32.** Torch trains in float32; an R
float64 forward on the exported (float32-valued) weights differs from
the torch embedding by \\\sim 10^{-6}\\. To guarantee the frozen linear
head is fit in EXACTLY the same embedding the scorer evaluates, the head
trains on the embeddings produced by the PURE-R forward over ALL
training rows (`.ssl_vicreg_forward`), NOT the torch embeddings. Fit and
score thus use the identical computation (the deepmaha / proto-net
fit/score-consistency lesson: fit and score must share one code path).

**The embedding forward.** The exported encoder is \\\mathrm{input}(p)
\to \[\mathrm{Linear} \to \phi\]\_{1} \to \dots \to \[\mathrm{Linear}
\to \phi\]\_{L-1} \to \mathrm{Linear} \to \mathrm{embedding}(d)\\. There
are `length(hidden)` linear layers up to the embedding (the last
`hidden` width is the embedding dim \\d\\); the activation (`"relu"`
default, or `"tanh"`) is applied AFTER every linear layer EXCEPT the
final embedding layer (no activation follows the embedding – it is the
raw representation, matching the torch encoder).
`.ssl_vicreg_forward(M, weights)` reimplements this exactly.

**Frozen linear-probe logit score.** With the frozen head weight \\w\\
(length \\d\\) and bias \\b\\, the score of an embedded specimen \\z\\
is the linear logit \$\$s(z) = w^\top z + b,\$\$ the SSL linear-probe
evaluation read-out. The head is a single BN-free \\\mathrm{Linear}(d
\to 1)\\ trained by seeded, fixed-epoch, `head_weight_decay`-regularised
optimisation on the FROZEN encoder embeddings (regularisation handles
separation on small \\n\\). Larger = more case-like.

**Single-sample transform.** Each specimen is mapped to the
self-contained per-sample robust CLR over its OWN strictly-positive
support (geometric-mean centring on `v > 0`) in the FROZEN feature
universe; this is exactly invariant to per-specimen scaling and uses no
cross-row statistic. Universe features absent from a specimen carry the
neutral rCLR value `0`. The SAME rCLR transform is applied to the frozen
training rows at fit and to every query at score, so a specimen's score
depends only on that specimen and the frozen model – it passes
[`singlesample_assert_row_equivariant`](https://kstawiski.github.io/OmicSelector/reference/singlesample_assert_row_equivariant.md).
The augmentation (masking + jitter) is applied ONLY at fit to form the
SSL views; scoring uses the clean rCLR with no augmentation and no RNG.

**Degenerate-neutral.** A query returns the neutral score `0` if fewer
than `hp$min_features` universe features overlap the query columns (a
column-overlap floor, batch-independent), or it has empty positive
support over the frozen universe tested on the (pre-rCLR) aligned
ORIGINAL abundances (`!any(X_use[i, ] > 0)`). A FLAT all-equal-positive
composition maps to the rCLR origin but is a VALID specimen and is
scored normally (NOT floored) – the empty-support test is on the
ORIGINAL abundances, NOT on `all(z == 0)`.

## References

Bardes A, Ponce J, LeCun Y. (2022) VICReg:
Variance-Invariance-Covariance Regularization for Self-Supervised
Learning. *International Conference on Learning Representations (ICLR)*.
arXiv:2105.04906.
