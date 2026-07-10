# Prototypical-network single-sample discriminator (frozen class prototypes)

Single-sample within-cohort discriminator (family N) implementing a
prototypical network in the sense of Snell, Swersky & Zemel (2017),
adapted to compositional specimens. A small BatchNorm-free MLP embedding
of the per-sample-rCLR specimen is learned on the TRAINING set by
EPISODIC PROTOTYPICAL training (NOT cross-entropy) via reticulate-python
torch, then **FROZEN**; per-class PROTOTYPES (class-mean embeddings) are
computed in the embedding space, and a new specimen is scored by the
difference of its squared-Euclidean distances to the frozen control and
case prototypes. Larger score = more case-like.

**Episodic prototypical training (Snell 2017).** Unlike the
cross-entropy-trained `lrt-deepmaha` sibling, the embedding here is
trained with the prototypical loss. Each episode draws a disjoint
SUPPORT set (`n_support` per class) and QUERY set (`n_query` per class)
from the (rCLR'd, universe-aligned) training rows; the prototype of
class \\c\\ is the mean SUPPORT embedding of that class; each query is
assigned \\\log\\\mathrm{softmax}(-d_c)\\ over its squared-Euclidean
distances \\d_c = \lVert z_q - \mathrm{proto}\_c\rVert^2\\ to the
prototypes, and the loss is the mean negative log-probability of the
query's TRUE class. This is exactly Snell 2017 eq. (1)-(2) with the
squared-Euclidean distance (the Bregman-divergence choice the paper
shows is principled). Adam optimises the embedding net over `epochs`
episodes.

**Torch is confined to FIT.** The MLP is trained with torch on the
rCLR'd, universe-aligned training matrix. After training, the weight
matrices and biases up to and INCLUDING the embedding layer are EXPORTED
to R as plain numeric matrices/vectors; the python module is DISCARDED.
The fitted model therefore holds NO external pointer and NO python
dependency at score time: the embedding forward, the frozen prototypes,
and the Euclidean-to-prototype score are all pure base-R. This makes the
default `model_digest` a stable deterministic snapshot and lets the
package suite run the scorer on any host without the venv.

**Fit/score consistency despite float32.** Torch trains in float32; an R
float64 forward on the exported (float32-valued) weights differs from
the torch embedding by \\\sim 10^{-6}\\. To guarantee that the frozen
prototypes live in EXACTLY the same embedding the scorer evaluates, the
prototypes \\\mu\_{\mathrm{case}}, \mu\_{\mathrm{ctrl}}\\ are computed
from the PURE-R forward over ALL training rows (`.proto_net_forward`),
NOT from the torch (or the episodic support-mean) embeddings. Fit and
score thus use the identical computation (the deepmaha / nc-ecod-copod /
nc-sinkhorn-single lesson: fit and score must share one code path).

**The embedding forward.** The exported net is \\\mathrm{input}(p) \to
\[\mathrm{Linear} \to \phi\]\_{1} \to \dots \to \[\mathrm{Linear} \to
\phi\]\_{L-1} \to \mathrm{Linear} \to \mathrm{embedding}(d)\\. There are
`length(hidden)` linear layers up to the embedding (the last `hidden`
width is the embedding dim \\d\\); the activation (`"relu"` default, or
`"tanh"`) is applied AFTER every linear layer EXCEPT the final embedding
layer (no activation follows the embedding – it is the raw
representation, matching the torch net).
`.proto_net_forward(M, weights)` reimplements this exactly.

**Euclidean-to-prototype score.** With the frozen prototypes the score
of an embedded specimen \\z\\ is \$\$s(z) = -\lVert z -
\mu\_{\mathrm{case}}\rVert^2 + \lVert z -
\mu\_{\mathrm{ctrl}}\rVert^2,\$\$ the squared-Euclidean log-likelihood
ratio to the prototypes. This is exactly the deepmaha score with
\\\Sigma = I\\ (no covariance / no inversion) – proto-net is the
identity-metric, prototype-mean special case, but with an embedding
trained by the PROTOTYPICAL loss, which is what makes it a distinct
method ("Euclidean-to-prototype metric learning"). Larger = more
case-like.

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
ORIGINAL abundances (`!any(X_use[i, ] > 0)`). A FLAT all-equal-positive
composition maps to the rCLR origin but is a VALID specimen and is
scored normally (NOT floored) – the empty-support test is on the
ORIGINAL abundances, NOT on `all(z == 0)`.

## References

Snell J, Swersky K, Zemel R. (2017) Prototypical Networks for Few-shot
Learning. *Advances in Neural Information Processing Systems (NeurIPS)*
30. arXiv:1703.05175.
