# DANN (domain-adversarial) cross-cohort transfer discriminator

Single-sample TRANSFER-estimand discriminator (family N) implementing
**DANN** (Domain-Adversarial Neural Network; Ganin et al. 2016, JMLR 17;
arXiv:1505.07818) for cross-cohort out-of-distribution generalization. A
small BatchNorm-free MLP encoder of the per-sample-rCLR specimen and a
linear classification (LABEL) head \\\mathrm{Linear}(d \to 1)\\ are
trained JOINTLY across multiple training ENVIRONMENTS (cohorts) by the
DANN objective via reticulate-python torch, then **FROZEN**; a new
specimen is scored by its frozen logit on its pure-R encoder embedding.
Larger score = more case-like.

**Transfer is in the FIT, not the score.** The cross-cohort aspect of
this method lives ENTIRELY in training: the encoder + label head are
optimised so the label loss is low while a DOMAIN discriminator, fed the
encoder embedding through a GRADIENT-REVERSAL LAYER (GRL), is unable to
recover the cohort. The reversed gradient drives the encoder toward a
representation that is INVARIANT to the cohort while remaining
label-predictive – DANN's mechanism for distributional robustness across
cohorts. The SCORE path is identical to the within-estimand
trained-torch siblings (`lrt-deepmaha`, `proto-net`, `ssl-vicreg`,
`dro-vrex`, `dg-fishr`): a frozen linear logit on a pure-R encoder
embedding of the single specimen, so single-row scoring equals batch
scoring EXACTLY and the score is exactly invariant to per-specimen
positive scaling.

**PRESPEC DEPARTURE (kit-conditional -\> cohort-adversarial).** The
roster names `dann` with the present kit-conditional
`os_fit_dann_kit_extended` (a kit/technology-EXTENDED design that
conditions the network on kit metadata). Kit metadata is an EXTERNAL
cross-sample covariate: feeding it at inference violates single-sample
deployability (no external covariate at score). This implements the
CANONICAL single-sample DANN instead: the encoder is trained to be
cohort-INVARIANT via a gradient-reversal domain discriminator at FIT, so
at SCORE only the frozen encoder + label head are used (NO kit/cohort
input). This is the standard DANN and is single-sample-faithful.

**The DANN objective (Ganin 2016).** Let \\z = \mathrm{encoder}(x)\\ be
the embedding. Three modules are trained JOINTLY: the encoder \\\phi\\,
a LABEL head \\\mathrm{Linear}(d \to 1)\\ with logit \\l = w^\top z +
b\\, and a DOMAIN discriminator \\D\\ (a small MLP \\d \to \dots \to K\\
over the \\K\\ cohorts). The total loss is \$\$\mathcal{L} =
\underbrace{\overline{\mathrm{BCE}}(l, y)}\_{\text{label loss}} \\+\\
\underbrace{\mathrm{CE}\big(D(\mathrm{GRL}\_\lambda(z)),\\
c\big)}\_{\text{domain loss}},\$\$ where \\c\\ is the cohort label and
\\\mathrm{GRL}\_\lambda\\ is the gradient-reversal layer: the IDENTITY
on the forward pass, and on the backward pass it multiplies the incoming
gradient by \\-\lambda\\ (\\\lambda = \\`grl_lambda`, default `1.0`).
Because the domain branch's gradient into the encoder is negated,
minimising \\\mathcal{L}\\ drives the discriminator to classify the
cohort WHILE the encoder is pushed to DEFEAT it – i.e. toward a
cohort-invariant \\z\\. The label head and encoder also minimise the
ordinary label loss, so \\z\\ stays label-predictive. Encoder + label
head + domain head are optimised jointly with Adam over `epochs`
full-batch steps.

**The GRL is FIT-only.** The gradient-reversal layer and the domain
discriminator exist ONLY during training (to shape \\\phi\\). They are
NOT exported and NOT used at score. The exported model holds ONLY the
encoder weights (up to and including the embedding) and the LABEL head
\\(w, b)\\.

**Environments from `meta_train`.** Environments (= domains for the
discriminator) are the distinct non-missing labels in
`meta_train[[cohort_col]]`. `cohort_col` defaults to `NULL` and is
auto-detected from common names (`"accession"`, `"cohort"`, `"batch"`,
`"study"`, `"dataset"`, `"group"`); set it explicitly to override. An
environment is VALID only if it carries at least one case AND one
control (a single-class environment has a degenerate,
label-uninformative risk and is DROPPED). The adversarial domain term
needs \\\ge 2\\ valid domains; if fewer remain (including `NULL`/absent
`meta_train`, a missing `cohort_col`, or a single distinct cohort),
fitting FALLS BACK to plain ERM (a single pooled label risk over all
rows, the domain term INERT, `grl_lambda` unused) – mirroring
`dro-vrex`'s and `dg-fishr`'s NULL-meta single-group ERM degeneration.
The realised number of domains is recorded as `model$n_environments`;
`n_environments == 1` means the fit ran as ERM and the adversary was NOT
engaged.

**Torch is confined to FIT.** After training, the encoder weight
matrices and biases up to and INCLUDING the embedding layer, and the
LABEL head \\(w, b)\\, are EXPORTED to R as plain numeric
matrices/vectors; the python module (encoder, label head, domain head,
GRL) and the torch objects are DISCARDED. The fitted model holds NO
external pointer and NO python dependency at score time: the encoder
forward and the linear-head logit are pure base-R. This makes the
default `model_digest` a stable deterministic snapshot and lets the
package suite run the scorer on any host without the venv.

**The embedding forward.** The exported encoder is \\\mathrm{input}(p)
\to \[\mathrm{Linear} \to \phi\]\_{1} \to \dots \to \[\mathrm{Linear}
\to \phi\]\_{L-1} \to \mathrm{Linear} \to \mathrm{embedding}(d)\\. There
are `length(hidden)` linear layers up to the embedding (the last
`hidden` width is the embedding dim \\d\\); the activation (`"relu"`
default, or `"tanh"`) is applied AFTER every linear layer EXCEPT the
final embedding layer. `.dann_forward(M, weights, activation)`
reimplements this exactly. The score of an embedded specimen \\z\\ is
the frozen linear LABEL logit \\s(z) = w^\top z + b\\; larger = more
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

Ganin Y, Ustinova E, Ajakan H, Germain P, Larochelle H, Laviolette F,
Marchand M, Lempitsky V. (2016) Domain-Adversarial Training of Neural
Networks. *Journal of Machine Learning Research* 17(59):1–35.
arXiv:1505.07818.
