# Domain-adversarial conditional VAE single-sample discriminator (DA-cVAE)

Single-sample within-cohort discriminator (family N) implementing Stage
2/3 of idea_report_1570797495812882433 ("Structurally Informed ILR +
DA-cVAE"). The encoder ingests the **Stage-1 structural-ILR balances**
(the genomic-cluster / GC-tertile Sequential Binary Partition of
[`fit_ss_struct_ilr`](https://kstawiski.github.io/OmicSelector/reference/fit_ss_struct_ilr.md),
computed per specimen via
[`ws_balance_ilr`](https://kstawiski.github.io/OmicSelector/reference/ws_balance_ilr.md)
and standardised with train-only statistics) and maps them to a latent
\\z\\. Training combines three FIT-only regularisers on that latent — a
conditional VAE decoder \\p(\mathrm{balances}\mid z,
\mathrm{platform})\\, a KL prior, and a **gradient-reversal platform
adversary** — with a supervised disease LABEL head. The encoder + label
head are then **FROZEN** and exported to R; a new specimen is scored in
PURE base-R by its frozen disease logit on the deterministic latent \\z
= \mu\\. Larger = more case-like.

**Why this architecture (cross-platform transportability).** The
platform/assay label conditions ONLY the decoder (conditional
generation) and drives the GRL adversary; the ENCODER never sees the
platform, so a held-out specimen is scored with no platform/cohort
covariate (single-sample-faithful). The reversed adversary gradient
pushes \\z\\ toward a representation whose platform is not recoverable —
the mechanism for *reduced platform imprint / improved cross-platform
transportability*. As enforced in the design (v2.1 §0), this is the
claim ceiling: the method targets reduced platform imprint, NOT platform
"invariance", and every result from it is EXPLORATORY and never merged
into the frozen primary headline.

**The DA-cVAE objective.** Let \\b\\ be the standardised structural-ILR
balances of a specimen, \\(\mu, \log\sigma^2) = \mathrm{enc}(b)\\, \\z =
\mu + \sigma \odot \epsilon\\ (\\\epsilon \sim N(0,I)\\; ONE draw per
row per step), \\e_g = \mathrm{onehot}(\mathrm{platform})\\. Four terms
are minimised jointly with Adam: \$\$\mathcal{L} =
\underbrace{\overline{\mathrm{BCE}}(w^\top\mu + b_0,\\ y)}\_{\text{label
(on }\mu\text{)}} + \rho\\\underbrace{\overline{\lVert \mathrm{dec}(z,
e_g) - b\rVert^2}}\_{\text{recon}} +
\beta\\\underbrace{\overline{\mathrm{KL}(q(z\mid
b)\\\Vert\\N(0,I))}}\_{\text{prior}} +
\underbrace{\mathrm{CE}\big(\mathrm{adv}(\mathrm{GRL}\_\lambda(z)),\\
g\big)}\_{\text{platform adversary}}.\$\$ The label head is trained on
the DETERMINISTIC \\\mu\\ (exactly the latent used at score), so the
kept readout sees no train/score representation mismatch. The decoder,
logvar head, and adversary are FIT-only and are NOT exported.

**Platforms from `meta_train`.** The platform/assay domains are the
distinct non-missing labels in `meta_train[[platform_col]]`
(auto-detected from `"platform"`, `"assay"`, `"technology"`,
`"accession"`, `"batch"`, `"study"`, `"dataset"` when `NULL`). The
conditional decoder uses ALL \\\ge 2\\ platforms whenever they exist,
INDEPENDENT of `grl_lambda`; the GRL adversary is engaged SEPARATELY,
only when `grl_lambda > 0` (recorded as `model$adv_engaged`). Setting
`grl_lambda = 0` with \\\ge 2\\ platforms is therefore the design's
**no-GRL ablation**: the SAME architecture (conditional decoder
preserved), adversary off — exactly the baseline the design's
non-negotiable no-harm guard compares against. A single platform (or
`NULL` `meta_train`) collapses the decoder to unconditional
(`n_platforms == 1`) with the adversary necessarily off.

**Torch is confined to FIT.** After training, the encoder weights (trunk
up to and including the \\\mu\\ head) and the disease LABEL head \\(w,
b_0)\\ are EXPORTED to R as plain float64 matrices/vectors; the python
module (logvar head, decoder, adversary, GRL) is DISCARDED. The fitted
model holds NO external pointer and NO python dependency at score time —
the encoder forward and the linear logit are pure base-R, so the suite
runs the scorer on any host without the venv, and the default
`model_digest` is a stable deterministic snapshot.

**The encoder forward + score.** The exported encoder is \\b \to
\[\mathrm{Linear}\to\phi\]\_{1}\to\dots\to\[\mathrm{Linear}\to\phi\]\_{H}\to
\mathrm{Linear}\to\mu(L)\\ (activation after every layer EXCEPT the
final \\\mu\\ layer; `"relu"` default or `"tanh"`), reimplemented
exactly by `.dacvae_forward`. The disease score is the frozen logit \\s
= w^\top\mu + b_0\\.

**Novelty (secondary capability).** A frozen kNN-distance novelty score
over the TRAIN control (reference) latents is also serialised (design
F9: multi-reference, not a single centroid);
[`score_dacvae_novelty`](https://kstawiski.github.io/OmicSelector/reference/score_dacvae_novelty.md)
returns, per query, the mean distance to its `novelty_k` nearest
reference latents. This supports the cross-platform /
out-of-distribution capability and is NOT the disease score.

**Single-sample contract.** Each specimen's balances are computed from
its OWN values (the
[`ws_balance_ilr`](https://kstawiski.github.io/OmicSelector/reference/ws_balance_ilr.md)
per-sample pseudocount makes them exactly invariant to per-specimen
positive scaling), standardised with the FROZEN train mean/sd (missing
balances imputed to the standardised mean `0`), and embedded by the
frozen encoder. No cross-row / batch statistic is used, so a specimen's
score depends only on that specimen and the frozen model — it passes
[`singlesample_assert_row_equivariant`](https://kstawiski.github.io/OmicSelector/reference/singlesample_assert_row_equivariant.md)
(single row == batch EXACTLY, exact per-specimen scale invariance).
Specimens with empty positive support over the frozen universe, or fewer
than `hp$min_features` universe features present, return the neutral
score `0`.

## References

Kingma DP, Welling M. (2014) Auto-Encoding Variational Bayes. *ICLR*.
arXiv:1312.6114.

Sohn K, Lee H, Yan X. (2015) Learning Structured Output Representation
using Deep Conditional Generative Models. *NeurIPS* 28.

Ganin Y, et al. (2016) Domain-Adversarial Training of Neural Networks.
*JMLR* 17(59):1–35. arXiv:1505.07818.
