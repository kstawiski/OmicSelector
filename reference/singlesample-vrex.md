# V-REx (Variance Risk Extrapolation) cross-cohort transfer discriminator

Single-sample TRANSFER-estimand discriminator (family N) implementing
**V-REx** (Variance Risk Extrapolation; Krueger et al. 2021 ICML,
arXiv:2003.00688) for cross-cohort out-of-distribution generalization. A
small BatchNorm-free MLP encoder of the per-sample-rCLR specimen and a
linear classification head \\\mathrm{Linear}(d \to 1)\\ are trained
JOINTLY across multiple training ENVIRONMENTS (cohorts) by the V-REx
objective via reticulate-python torch, then **FROZEN**; a new specimen
is scored by its frozen logit on its pure-R encoder embedding. Larger
score = more case-like.

**Transfer is in the FIT, not the score.** The cross-cohort aspect of
this method lives ENTIRELY in training: the encoder + head are optimised
so the per-environment risks are both LOW (mean term) and SIMILAR
(variance term), which is V-REx's mechanism for distributional
robustness across cohorts. The SCORE path is identical to the
within-estimand trained-torch siblings (`lrt-deepmaha`, `proto-net`,
`ssl-vicreg`): a frozen linear logit on a pure-R encoder embedding of
the single specimen, so single-row scoring equals batch scoring EXACTLY
and the score is exactly invariant to per-specimen positive scaling.

**The V-REx objective (Krueger 2021).** For environments \\e \in
\mathcal{E}\\, let \\R_e\\ be the mean binary-cross-entropy-with-logits
risk on environment \\e\\'s rows (label vs head logit on the encoder
embedding). V-REx minimises \$\$\mathcal{L} =
\frac{1}{\|\mathcal{E}\|}\sum_e R_e \\+\\ \beta \cdot
\mathrm{Var}\_e(R_e),\$\$ where \\\mathrm{Var}\_e(R_e) =
\frac{1}{\|\mathcal{E}\|}\sum_e (R_e - \bar R)^2\\ is the BIASED
(population, \\1/\|\mathcal{E}\|\\) variance of the per-environment
risks and \\\beta = \\`vrex_beta` (default `1.0`). Encoder + head are
optimised jointly with Adam over `epochs` full-batch steps. The variance
penalty equalises risk across cohorts: a predictor that is accurate in
one cohort but not another is penalised, biasing the solution toward
features that generalise.

**Environments from `meta_train`.** Environments are the distinct
non-missing labels in `meta_train[[cohort_col]]`. `cohort_col` defaults
to `NULL` and is auto-detected from common names (`"accession"`,
`"cohort"`, `"batch"`, `"study"`, `"dataset"`, `"group"`); set it
explicitly to override. An environment is VALID for the V-REx penalty
only if it carries at least one case AND one control (a single-class
environment has a degenerate, label-uninformative risk and is DROPPED
from the per-environment risk set). If fewer than two valid environments
remain (including `NULL`/absent `meta_train`, a missing `cohort_col`, or
a single distinct cohort), fitting FALLS BACK to plain ERM (a single
pooled risk over all rows, `vrex_beta` inert) – mirroring `dro-group`'s
NULL-meta single-group ERM degeneration. The realised number of
environments is recorded as `model$n_environments`;
`n_environments == 1` means the fit ran as ERM and the transfer estimand
was NOT engaged.

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
final embedding layer. `.dro_vrex_forward(M, weights, activation)`
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

Krueger D, Caballero E, Jacobsen J-H, Zhang A, Binas J, Zhang D, Le
Priol R, Courville A. (2021) Out-of-Distribution Generalization via Risk
Extrapolation (REx). *Proceedings of the 38th International Conference
on Machine Learning (ICML)*, PMLR 139. arXiv:2003.00688.
