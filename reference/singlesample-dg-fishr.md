# Fishr (gradient-variance matching) cross-cohort transfer discriminator

Single-sample TRANSFER-estimand discriminator (family N) implementing
**Fishr** (Rame, Dancette, Cord 2022 ICML; arXiv:2109.02934) for
cross-cohort out-of-distribution generalization. A small BatchNorm-free
MLP encoder of the per-sample-rCLR specimen and a linear classification
head \\\mathrm{Linear}(d \to 1)\\ are trained JOINTLY across multiple
training ENVIRONMENTS (cohorts) by the Fishr objective via
reticulate-python torch, then **FROZEN**; a new specimen is scored by
its frozen logit on its pure-R encoder embedding. Larger score = more
case-like.

**Transfer is in the FIT, not the score.** The cross-cohort aspect of
this method lives ENTIRELY in training: the encoder + head are optimised
so the per-environment risks are low (mean term) and the per-environment
classifier-gradient VARIANCES are matched (Fishr term), which is Fishr's
mechanism for distributional robustness across cohorts. The SCORE path
is identical to the within-estimand trained-torch siblings
(`lrt-deepmaha`, `proto-net`, `ssl-vicreg`, `dro-vrex`): a frozen linear
logit on a pure-R encoder embedding of the single specimen, so
single-row scoring equals batch scoring EXACTLY and the score is exactly
invariant to per-specimen positive scaling.

**The Fishr objective (Rame 2022).** Let \\z = \mathrm{encoder}(x)\\ be
the embedding and the head logit be \\l = w^\top z + b\\. The per-sample
gradient of the binary-cross-entropy-with-logits loss w.r.t. the head
pre-activation logit is the classic residual \\g_i = \sigma(l_i) -
y_i\\; the per-sample gradient w.r.t. the head WEIGHTS is \\G_i = g_i \\
z_i\\ (a length-\\d\\ vector; the bias gradient \\g_i\\ is appended,
giving a length-\\d{+}1\\ per-sample gradient vector). For each VALID
environment \\e\\, the per-environment gradient-variance vector is the
diagonal of the per-env gradient covariance, \$\$v_e = \mathrm{Var}\_{i
\in e}(G_i) = \frac{1}{n_e}\sum\_{i \in e} (G_i - \bar G_e)^2
\quad(\text{biased, } 1/n_e),\$\$ and the Fishr penalty is the
across-environment variance of those vectors, \$\$\Omega =
\frac{1}{\|\mathcal{E}\|}\sum_e \lVert v_e - \bar v \rVert_2^2, \qquad
\bar v = \frac{1}{\|\mathcal{E}\|}\sum_e v_e.\$\$ Fishr minimises
\\\mathcal{L} = \overline{\mathrm{BCE}} + \lambda\\\Omega\\ with
\\\lambda = \\`fishr_lambda` (default `1.0`);
\\\overline{\mathrm{BCE}}\\ is the mean per-row risk over the retained
valid-environment rows. Encoder + head are optimised jointly with Adam
over `epochs` full-batch steps. Matching the per-environment gradient
variances drives the classifier toward features whose loss-landscape
geometry is consistent across cohorts.

**Environments from `meta_train`.** Environments are the distinct
non-missing labels in `meta_train[[cohort_col]]`. `cohort_col` defaults
to `NULL` and is auto-detected from common names (`"accession"`,
`"cohort"`, `"batch"`, `"study"`, `"dataset"`, `"group"`); set it
explicitly to override. An environment is VALID for the Fishr penalty
only if it carries at least one case AND one control (a single-class
environment has a degenerate, label-uninformative risk and is DROPPED
from the per-environment gradient set). If fewer than two valid
environments remain (including `NULL`/absent `meta_train`, a missing
`cohort_col`, or a single distinct cohort), fitting FALLS BACK to plain
ERM (a single pooled risk over all rows, `fishr_lambda` inert) –
mirroring `dro-vrex`'s and `dro-group`'s NULL-meta single-group ERM
degeneration. The realised number of environments is recorded as
`model$n_environments`; `n_environments == 1` means the fit ran as ERM
and the transfer estimand was NOT engaged.

**Torch is confined to FIT.** After training, the encoder weight
matrices and biases up to and INCLUDING the embedding layer, and the
head \\(w, b)\\, are EXPORTED to R as plain numeric matrices/vectors;
the python module and the torch objects are DISCARDED. The fitted model
holds NO external pointer and NO python dependency at score time: the
encoder forward and the linear-head logit are pure base-R. This makes
the default `model_digest` a stable deterministic snapshot and lets the
package suite run the scorer on any host without the venv.

**The embedding forward.** The exported encoder is \\\mathrm{input}(p)
\to \[\mathrm{Linear} \to \phi\]\_{1} \to \dots \to \[\mathrm{Linear}
\to \phi\]\_{L-1} \to \mathrm{Linear} \to \mathrm{embedding}(d)\\. There
are `length(hidden)` linear layers up to the embedding (the last
`hidden` width is the embedding dim \\d\\); the activation (`"relu"`
default, or `"tanh"`) is applied AFTER every linear layer EXCEPT the
final embedding layer. `.dg_fishr_forward(M, weights, activation)`
reimplements this exactly. The score of an embedded specimen \\z\\ is
the frozen linear logit \\s(z) = w^\top z + b\\; larger = more
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

Rame A, Dancette C, Cord M. (2022) Fishr: Invariant Gradient Variances
for Out-of-Distribution Generalization. *Proceedings of the 39th
International Conference on Machine Learning (ICML)*, PMLR 162.
arXiv:2109.02934.
