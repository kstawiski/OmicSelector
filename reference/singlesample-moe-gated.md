# Self-gated mixture of frozen experts (MoE-gated), single-sample

Single-sample within-cohort discriminator (family N) implementing a
**self-gated mixture of experts** over the specimen's OWN rCLR
embedding, adapted to compositional specimens. A small **LayerNorm** MLP
encoder of the per-sample-rCLR specimen is trained on the TRAINING set
(via reticulate-python torch) to a pre-activation (penultimate)
embedding \\z\\, then **FROZEN**; \\E\\ linear logistic expert heads
each give an expert logit \\h_e(z)\\, and a linear gating network with a
softmax over the \\E\\ experts gives mixture weights \\g_e(z)\\; the
score is the mixture logit \\s(z) = \sum_e g_e(z)\\h_e(z)\\. Larger
score = more case-like.

**Single-sample faithfulness: SELF-gated, NOT kit-gated.** The manifest
names this method "kit-gated mixture of frozen experts," where the gate
is driven by sequencing-kit / technology metadata. Kit/technology is an
EXTERNAL, cross-sample covariate; conditioning the gate on it at
inference would require an external label per query and thereby VIOLATE
single-sample deployability (no cross-cohort reference, no external
covariate at score). The single-sample- faithful design here is
therefore **self-gated**: the gate is a function of the specimen's OWN
frozen embedding \\z\\ alone, so a new specimen is routed purely from
its own profile with no external signal. This is the design departure
flagged for the gate; the MoE structure (experts + a learned soft
router) is preserved, only the routing SIGNAL is made endogenous.

**Torch is confined to FIT.** The encoder, the \\E\\ experts, and the
gate are trained JOINTLY with torch (Adam, cross-entropy on the mixture
logit) on the rCLR'd, universe-aligned training matrix. After training,
the encoder Linear weights/biases, the per-layer LayerNorm
\\\gamma,\beta,\varepsilon\\, the \\E\\ expert weight/bias and the gate
weight/bias are EXPORTED to R as plain numeric matrices/vectors; the
python module is DISCARDED. The fitted model therefore holds NO external
pointer and NO python dependency at score time: the encoder forward
(incl. LayerNorm), the expert logits, the gate softmax, and the mixture
are all pure base-R. This makes the default `model_digest` a stable
deterministic snapshot and lets the package suite run the scorer on any
host without the venv.

**Pre-activation embedding + LayerNorm (single-sample-safe).** The
encoder embeds to the PRE-ACTIVATION (penultimate) representation \\z\\
(the endorsed trained-torch lesson): activation is applied after every
LayerNorm EXCEPT the final (embedding) layer, so \\z\\ is the raw
representation the experts and the gate read. **LayerNorm** normalises
ACROSS FEATURES WITHIN A ROW, so it is per-row and single-sample-exact;
**BatchNorm**, which couples rows, is forbidden here. torch's LayerNorm
uses the BIASED variance (1/N), reproduced R-side with `rowMeans` of
squared deviations.

**Fit/score consistency despite float32.** Torch trains in float32; an R
float64 forward on the exported (float32-valued) weights differs from
the torch forward by \\\sim 10^{-6}\\. The exported experts and gate are
torch-trained jointly with the encoder, so the mixture is read out from
the SAME exported weights the scorer evaluates; the pure-R forward,
expert logits, and gate softmax are the deployed computation (the
deepmaha / proto-net lesson: fit and score share one code path).

**The mixture forward.** With the frozen embedding \\z =
\mathrm{enc}(\mathrm{rclr}(x))\\, the \\E\\ expert logits are \\H =
z\\W_e^\top + b_e\\ (length \\E\\), the gate logits are \\z\\W_g^\top +
b_g\\, the mixture weights are the softmax \\g_e(z) =
\mathrm{softmax}(z\\W_g^\top + b_g)\_e\\ (a valid distribution over the
experts, rows sum to 1, all \\\geq 0\\), and the score is the convex
combination \\s(z) = \sum_e g_e(z)\\h_e(z)\\. Because the gate is a
softmax, the mixture logit lies within \\\[\min_e h_e, \max_e h_e\]\\
per row. Everything is per-row, so single-row == batch EXACTLY.

**Single-sample transform.** Each specimen is mapped to the
self-contained per-sample robust CLR over its OWN strictly-positive
support (geometric-mean centring on `v > 0`) in the FROZEN feature
universe; this is exactly invariant to per-specimen scaling and uses no
cross-row statistic. Universe features absent from a specimen carry the
neutral rCLR value `0`. The SAME rCLR transform is applied to the frozen
training rows at fit and to every query at score; combined with
LayerNorm (per-row), the frozen experts, and the frozen gate, a
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

Jacobs RA, Jordan MI, Nowlan SJ, Hinton GE. (1991) Adaptive Mixtures of
Local Experts. *Neural Computation* 3(1):79-87. Shazeer N, Mirhoseini A,
Maziarz K, Davis A, Le Q, Hinton G, Dean J. (2017) Outrageously Large
Neural Networks: The Sparsely-Gated Mixture-of-Experts Layer.
*International Conference on Learning Representations (ICLR)*.
arXiv:1701.06538.
