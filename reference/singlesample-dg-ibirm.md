# IB-IRM (Information-Bottleneck Invariant Risk Minimization) cross-cohort transfer discriminator

Single-sample TRANSFER-estimand discriminator (family N) implementing
**IB-IRM** (Ahuja, Caballero, Zhang et al. 2021 NeurIPS;
arXiv:2106.06607) for cross-cohort out-of-distribution generalization. A
small BatchNorm-free MLP encoder of the per-sample-rCLR specimen and a
linear classification head \\\mathrm{Linear}(d \to 1)\\ are trained
JOINTLY across multiple training ENVIRONMENTS (cohorts) by the IB-IRM
objective via reticulate-python torch, then **FROZEN**; a new specimen
is scored by its frozen logit on its pure-R encoder embedding. Larger
score = more case-like.

**Transfer is in the FIT, not the score.** The cross-cohort aspect of
this method lives ENTIRELY in training: the encoder + head are optimised
so the per-environment risks are low (mean term), the per-environment
classifier is SIMULTANEOUSLY OPTIMAL across cohorts (the IRMv1
gradient-penalty term), and the representation is COMPRESSED (the
Information-Bottleneck embedding-variance term), which is IB-IRM's
mechanism for distributional robustness across cohorts. The SCORE path
is identical to the within-estimand trained-torch siblings
(`lrt-deepmaha`, `proto-net`, `ssl-vicreg`, `dro-vrex`, `dg-fishr`): a
frozen linear logit on a pure-R encoder embedding of the single
specimen, so single-row scoring equals batch scoring EXACTLY and the
score is exactly invariant to per-specimen positive scaling.

**The IB-IRM objective (Ahuja 2021).** Let \\z = \mathrm{encoder}(x)\\
be the embedding and the head logit be \\l = w_0^\top z + b_0\\. IB-IRM
minimises \$\$\mathcal{L} = \frac{1}{\|\mathcal{E}\|}\sum_e R_e \\+\\
\lambda\_{\mathrm{IRM}} \cdot \Omega\_{\mathrm{IRM}} \\+\\
\lambda\_{\mathrm{IB}} \cdot \Omega\_{\mathrm{IB}},\$\$ where \\R_e\\ is
the mean binary-cross-entropy-with-logits risk on environment \\e\\. The
two penalties are:

1.  **IRMv1 gradient penalty** (Arjovsky 2019, the IRM term IB-IRM
    inherits). A FIXED dummy scalar classifier multiplier \\w = 1.0\\ is
    applied to the head logit (\\l_e = w \cdot (w_0^\top z_e + b_0)\\);
    for each VALID environment \\e\\, \\R_e(w)\\ is the mean
    BCE-with-logits at \\w\\, and \$\$\Omega\_{\mathrm{IRM}} =
    \frac{1}{\|\mathcal{E}\|}\sum_e \left( \left. \frac{\partial
    R_e}{\partial w} \right\|\_{w=1} \right)^2.\$\$ The gradient
    \\\partial R_e / \partial w\\ is taken by torch
    `autograd.grad(..., create_graph = TRUE)` so the penalty itself
    stays differentiable w.r.t. the encoder + head. Because \\w\\ is a
    scalar, \\\lVert \partial R_e / \partial w \rVert^2\\ is just the
    SQUARE of that scalar gradient. This is the canonical IRMv1 penalty.

2.  **Information-Bottleneck (IB) variance penalty** (the term that
    distinguishes IB-IRM from plain IRM). The representation is
    compressed by penalising the per-dimension variance of the embedding
    over the POOLED retained rows: \$\$\Omega\_{\mathrm{IB}} =
    \frac{1}{d}\sum\_{j=1}^{d} \mathrm{Var}\_i(z\_{ij}) =
    \frac{1}{d}\sum\_{j=1}^{d} \frac{1}{n}\sum\_{i=1}^{n}\left(z\_{ij} -
    \bar z\_{\cdot j}\right)^2,\$\$ the mean over embedding dimensions
    of the BIASED (population, \\1/n\\) variance of each embedding
    feature. Minimising representation variance is the standard
    tractable IB surrogate (lower feature variance \\\to\\ lower
    representation entropy \\\to\\ a more compressed, invariant code).
    The biased \\1/n\\ variance is computed EXPLICITLY (not `torch$var`,
    whose default is the \\N-1\\ unbiased estimator).

Both \\\lambda\_{\mathrm{IRM}}\\ (`irm_lambda`, default `1.0`) and
\\\lambda\_{\mathrm{IB}}\\ (`ib_lambda`, default `1.0`) are
non-negative. Encoder + head are optimised jointly with Adam over
`epochs` full-batch steps.

**Environments from `meta_train`.** Environments are the distinct
non-missing labels in `meta_train[[cohort_col]]`. `cohort_col` defaults
to `NULL` and is auto-detected from common names (`"accession"`,
`"cohort"`, `"batch"`, `"study"`, `"dataset"`, `"group"`); set it
explicitly to override. An environment is VALID for the IRM penalty only
if it carries at least one case AND one control (a single-class
environment has a degenerate, label-uninformative risk and is DROPPED
from the per-environment gradient set). If fewer than two valid
environments remain (including `NULL`/absent `meta_train`, a missing
`cohort_col`, or a single distinct cohort), fitting FALLS BACK to plain
ERM (a single pooled risk over all rows, BOTH penalties inert) –
mirroring `dro-vrex`'s and `dg-fishr`'s NULL-meta single-environment ERM
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
final embedding layer. `.dg_ibirm_forward(M, weights, activation)`
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

Ahuja K, Caballero E, Zhang D, Gagnon-Audet J-C, Bengio Y, Mitliagkas I,
Rish I. (2021) Invariance Principle Meets Information Bottleneck for
Out-of-Distribution Generalization. *Advances in Neural Information
Processing Systems (NeurIPS) 34*. arXiv:2106.06607.

Arjovsky M, Bottou L, Gulrajani I, Lopez-Paz D. (2019) Invariant Risk
Minimization. arXiv:1907.02893.
