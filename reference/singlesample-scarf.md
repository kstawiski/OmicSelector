# SCARF self-supervised contrastive single-sample discriminator (frozen linear head)

Single-sample within-cohort discriminator (family N) implementing the
SCARF self-supervised contrastive representation of Bahri, Jiang, Tay &
Metzler (2022), adapted to compositional specimens. A small
BatchNorm-free MLP encoder of the per-sample-rCLR specimen is learned on
the TRAINING set by SELF-SUPERVISED CONTRASTIVE training (InfoNCE /
NT-Xent, NO labels) via reticulate-python torch, then **FROZEN**; a
FROZEN linear-probe head (the SCARF evaluation protocol) is fit on the
labeled training embeddings, and a new specimen is scored by its frozen
linear-head logit on its pure-R encoder embedding. Larger score = more
case-like.

**SCARF contrastive training (Bahri 2022).** Unlike the cross-entropy
`lrt-deepmaha` or the prototypical `proto-net` siblings, the encoder
here is trained by the SCARF InfoNCE loss WITHOUT labels. For a batch of
the (rCLR'd, universe-aligned) training rows, a CORRUPTED view \\\tilde
x\\ of each row is built (SCARF eq. 1): independently for each feature,
with probability `corrupt_p` the entry is replaced by a value sampled
from that feature's EMPIRICAL TRAINING MARGINAL (sampled uniformly from
the training column's observed rCLR values). The positive pair is the
row and its corrupted view, both embedded through encoder then
projection head (\\z_i = g(f(x_i))\\, \\\tilde z_i = g(f(\tilde
x_i))\\); the NT-Xent contrastive loss with cosine similarity and
temperature \\\tau\\ has the positive pair \\(z_i, \tilde z_i)\\ as the
numerator and the other batch items as negatives, with the partition
function summing over ALL \\k\\ (including the positive \\k = i\\, per
SCARF Algorithm 1): \$\$\ell = \frac1B \sum_i -\log
\frac{\exp(\mathrm{sim}(z_i,\tilde z_i)/\tau)}
{\sum\_{k}\exp(\mathrm{sim}(z_i,\tilde z_k)/\tau)}.\$\$ This is exactly
softmax cross-entropy with diagonal labels. The corruption uses an
INDEPENDENT `Bernoulli(corrupt_p)` mask per feature (expected \\\lfloor
corrupt\\p \cdot M\rfloor\\ corrupted features per row), a common
variant of SCARF Algorithm 1's fixed-count \\q = \lfloor corrupt\\p
\cdot M\rfloor\\ per-sample selection. Adam optimises encoder +
projection head over `epochs` steps. The projection head is discarded
after training (standard contrastive practice); only the ENCODER (the
pre-projection representation) is exported.

**The corruption marginal is FIT-ONLY.** The per-feature empirical
marginal (the frozen training column) is used ONLY at FIT, to build
training views. It is NEVER stored in the model and NEVER consulted at
score: the score does a CLEAN encoder forward on the query, with no
corruption. Single-sample deployability therefore holds – a query's
score depends only on that query and the frozen encoder + frozen linear
head, with no cross-row statistic. Corruption randomness (the Bernoulli
mask and the marginal-resample row indices) is drawn via torch's SEEDED
generator, so two seed-matched CPU fits build bit-identical training
views.

**Frozen linear head (SCARF evaluation protocol).** After SSL training
the encoder is FROZEN. A single BN-free torch `Linear(d -> 1)` head is
then fit by seeded, fixed-epoch, `head_weight_decay`-regularised
binary-cross-entropy optimisation on the frozen encoder embeddings of
the labeled training rows (the L2 regularisation handles separation on
small \\n\\); its weight \\w\\ (length \\d\\) and bias \\b\\ (scalar)
are EXPORTED to R. CRITICAL fit/score consistency: the embeddings the
head trains on are computed from the SAME pure-R `.ai_scarf_forward` the
scorer uses (NOT the float32 torch embedding), so the head lives in
exactly the embedding the scorer evaluates (the deepmaha / proto-net
lesson: fit and score share one embedding code path).

**Torch is confined to FIT.** The encoder and head are trained with
torch on the rCLR'd, universe-aligned training matrix. After training,
the encoder weight matrices/biases and the linear-head \\w, b\\ are
EXPORTED to R as plain numeric matrices/vectors; the projection head and
the python module are DISCARDED. The fitted model therefore holds NO
external pointer and NO python dependency at score time: the encoder
forward and the linear-head logit are pure base-R. This makes the
default `model_digest` a stable deterministic snapshot and lets the
package suite run the scorer on any host without the venv.

**The encoder forward.** The exported net is \\\mathrm{input}(p) \to
\[\mathrm{Linear} \to \phi\]\_{1} \to \dots \to \[\mathrm{Linear} \to
\phi\]\_{L-1} \to \mathrm{Linear} \to \mathrm{embedding}(d)\\. There are
`length(hidden)` linear layers up to the embedding (the last `hidden`
width is the embedding dim \\d\\); the activation (`"relu"` default, or
`"tanh"`) is applied AFTER every linear layer EXCEPT the final embedding
layer (no activation follows the embedding – it is the raw
representation, matching the torch encoder).
`.ai_scarf_forward(M, weights)` reimplements this exactly.

**Linear-head logit score.** With the frozen head the score of an
embedded specimen \\z\\ is \\s(z) = w \cdot z + b\\, the linear-probe
logit. Larger = more case-like. Encoder forward + linear head are
per-row with no cross-row coupling, so the single-row score equals the
batch score EXACTLY.

**Single-sample transform.** Each specimen is mapped to the
self-contained per-sample robust CLR over its OWN strictly-positive
support (geometric-mean centring on `v > 0`) in the FROZEN feature
universe; this is exactly invariant to per-specimen scaling and uses no
cross-row statistic. Universe features absent from a specimen carry the
neutral rCLR value `0`. The SAME rCLR transform is applied to the frozen
training rows at fit and to every query at score, so a specimen's score
depends only on that specimen and the frozen model – it passes
[`singlesample_assert_row_equivariant`](https://kstawiski.github.io/OmicSelector/reference/singlesample_assert_row_equivariant.md).

**eval-frozen-BN.** The encoder is BN-FREE (the cleanest, exact
single-sample forward – same as deepmaha / proto-net). The manifest's
"eval-frozen-BN" alternative is acceptable only if BN running stats are
frozen and exported so the forward is per-row independent; BN-free is
simpler and exactly per-row, so this method uses BN-free.

**Degenerate-neutral.** A query returns the neutral score `0` if fewer
than `hp$min_features` universe features overlap the query columns (a
column-overlap floor, batch-independent), or it has empty positive
support over the frozen universe tested on the (pre-rCLR) aligned
ORIGINAL abundances (`!any(X_use[i, ] > 0)`). A FLAT all-equal-positive
composition maps to the rCLR origin but is a VALID specimen and is
scored normally (NOT floored) – the empty-support test is on the
ORIGINAL abundances, NOT on `all(z == 0)`.

## References

Bahri D, Jiang H, Tay Y, Metzler D. (2022) SCARF: Self-Supervised
Contrastive Learning using Random Feature Corruption. *International
Conference on Learning Representations (ICLR)*. arXiv:2106.15147.
