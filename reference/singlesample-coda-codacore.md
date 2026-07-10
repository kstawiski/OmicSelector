# CoDaCoRe stagewise log-ratio-balance single-sample discriminator (torch relaxation at fit, frozen discrete balances, PURE-R score)

Single-sample within-cohort discriminator (family O, compositional)
implementing CoDaCoRe (Gordon-Rodriguez, Quinn & Cunningham 2021,
"Learning sparse log-ratios for high-throughput sequencing data"): a
STAGEWISE (boosting) sequence of sparse log-ratio BALANCES, each learned
by a CONTINUOUS RELAXATION, scored through a FROZEN logistic head.
Larger score = more case-like.

**The CoDaCoRe algorithm.** For \\b = 1, \dots, B\\ (B =
`hp$max_balances`) a stage introduces a learnable per-feature assignment
over three states – numerator (+), denominator (-), and neither – via a
temperature-`tau` softmax of a logit matrix \\W_b \in \mathbb{R}^{p
\times 3}\\. The SOFT balance of specimen \\i\\ is the relaxed
log-contrast \$\$SB_b(i) = \frac{\sum_j s^{+}\_{j}\\\log x\_{ij}}{\sum_j
s^{+}\_{j}} - \frac{\sum_j s^{-}\_{j}\\\log x\_{ij}}{\sum_j
s^{-}\_{j}},\$\$ with \\s^{+}, s^{-}\\ the soft (+)/(-) memberships.
\\W_b\\ and a scalar stage slope are trained by Adam to minimise the
binary cross-entropy of the boosting RESIDUAL (\\y -
\sigma(\text{running score})\\). The trained soft memberships are then
DISCRETIZED to two disjoint feature sets \\A_b\\ (numerator, argmax
\\=+\\) and \\B_b\\ (denominator, argmax \\=-\\), with explicit floors
\\\|A_b\| \ge 1\\ and \\\|B_b\| \ge 1\\ and a deterministic tie-break;
\\A_b, B_b\\ are FROZEN. The discrete SBP/ILR balance is
\$\$\mathrm{bal}\_b(x) = \sqrt{\frac{\|A_b\|\|B_b\|}{\|A_b\|+\|B_b\|}}\\
\bigl(\overline{\log x\_{A_b}} - \overline{\log x\_{B_b}}\bigr),\$\$
computed over the specimen's OWN strictly-positive support within \\A_b
/ B_b\\. The running score is re-fit (a base-R IRLS ridge logistic on
the discrete-balance matrix), the residual recomputed, and the next
stage proceeds; a stage is ACCEPTED only if it lowers the training-only
deviance, otherwise stagewise growth stops.

**Distinct from selbal.** `bal-selbal` learns ONE balance by GREEDY
forward feature selection; `coda-codacore` learns a SEQUENCE of balances
by CONTINUOUS-RELAXATION stagewise boosting – the relaxation +
discretization + boosting residual are what make it CoDaCoRe, not
selbal.

**DEP-ROUTE PIVOT: codacore-tensorflow -\> reticulate-torch-plus-R.**
The prespecified backend was TensorFlow/Keras, but TensorFlow has NO
py3.14 wheel and the R codacore package is not installed. The CoDaCoRe
ALGORITHM is preserved EXACTLY; only the autodiff backend changes
TF/Keras -\> the already-installed torch (the relaxation is a tiny
length-\\3p\\ logit model). Torch is confined to FIT: the relaxation is
trained, the balances DISCRETIZED + FROZEN, and the discrete sets +
logistic head EXPORTED to R; the python module is DISCARDED. The fitted
model holds NO external pointer and NO python dependency at score time.

**The score is PURE BASE R.** For a query row \\x\\: align to the frozen
universe; for each frozen balance \\b\\ compute \\\mathrm{bal}\_b(x)\\
over \\x\\'s own strictly-positive support in \\A_b / B_b\\; if a side
has NO positive feature for this specimen that balance contributes its
neutral value `0`. The score is \\\mathrm{intercept} + \sum_b
w_b\\\mathrm{bal}\_b(x)\\. Each balance is a per-specimen difference of
mean-logs over that specimen's own features – it uses no cross-sample
statistic and no python – so a row's score is bit-identical scored alone
or in any batch (single-row == batch EXACTLY 0), and is exactly
invariant to per-sample positive scaling (a positive rescale \\c\\
cancels: \\\overline{\log(c\\x_A)} - \overline{\log(c\\x_B)} =
\overline{\log x_A} - \overline{\log x_B}\\).

**Reproducibility: the relaxation is fit on CPU by default.** The
DISCRETE balances (the \\A_b, B_b\\ sets) must be identical across two
seed-matched fits for the determinism gate and the `model_digest`. GPU
reductions are nondeterministic and could flip a near-threshold
feature's argmax membership, so the relaxation is trained on CPU by
DEFAULT (the per-stage model is tiny, so CPU is fast and EXACTLY
reproducible), seeded before build, full-batch with no shuffle, with a
deterministic argmax discretization and an explicit lowest-index
tie-break. (Tests may probe `device = "cuda"`/`"auto"` for the FIT; the
SCORE is pure-R and therefore device-irrelevant.)

**Degenerate-neutral / rCLR-origin.** A query returns the neutral score
`0` only if fewer than `hp$min_features` universe features overlap the
query columns (a column-overlap floor, batch-independent) OR it has
empty positive support over the frozen universe on the (pre-log) aligned
ORIGINAL abundances (`!any(X_use[i, ] > 0)`). A FLAT all-equal-positive
composition gives every balance the value `0` (\\\overline{\log x_A} -
\overline{\log x_B} = 0\\) but is a VALID specimen, scored to the
intercept (NOT floored). The score consumes NO RNG.

## References

Gordon-Rodriguez E, Quinn TP, Cunningham JP. (2021) Learning sparse
log-ratios for high-throughput sequencing data. *Bioinformatics*
38(1):157-163.

Egozcue JJ, Pawlowsky-Glahn V. (2005) Groups of parts and their balances
in compositional data analysis. *Mathematical Geology* 37(7):795-828.
