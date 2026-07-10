# DeepCoDA zero-sum log-contrast bottleneck + self-explaining head (single-sample)

Single-sample within-cohort discriminator (family O, compositional)
implementing DeepCoDA in the sense of Quinn, Nguyen, Rana, Gupta &
Venkatesh (2020), adapted to the single-specimen scoring contract. A
small set of ZERO-SUM LOG-CONTRAST bottlenecks (scale-invariant balances
/ log-ratios of the per-sample-rCLR specimen) and a SELF-EXPLAINING
per-sample predictor over those bottlenecks are learned end-to-end on
the TRAINING set by binary cross-entropy (via reticulate-python torch),
then **FROZEN**; a new specimen is scored by its frozen DeepCoDA logit.
Larger score = more case-like.

**Zero-sum log-contrast bottleneck (DeepCoDA eq. 1).** The first layer
maps the (per-sample rCLR'd, universe-aligned) input \\x\\ of length
\\p\\ to \\B\\ bottleneck "balances" \$\$b_k(x) = \sum_j W\_{kj}\\
\mathrm{rclr}(x)\_j,\$\$ under the HARD CONSTRAINT \\\sum_j W\_{kj} =
0\\ for every bottleneck \\k\\ (the log-contrast constraint – each
bottleneck is a valid log-ratio, invariant to the specimen total). The
constraint is enforced at training by parameterising \\W =
W^{\mathrm{raw}} - \mathrm{rowMeans}(W^{\mathrm{raw}})\\ (each row is
re-centred to sum to zero); the centred (zero-sum) \\W\\ is EXPORTED.
There is NO activation on the bottleneck (it is a linear log-contrast),
which also gives the pre-activation-embedding consistency of the
deepmaha / proto-net / vicreg siblings. Because a zero-sum contrast
annihilates any per-row additive constant, applying it to the rCLR is
identical to applying it to the raw log-abundance, so the bottleneck –
and therefore the whole score – is EXACTLY per-sample scale-invariant
(verified maxdiff \\\sim 0\\, up to float).

**Self-explaining per-sample predictor (Quinn 2020).** DeepCoDA's
"personalised interpretability" head is a self-explaining predictor: a
small BatchNorm-free MLP \\\theta(b(x))\\ produces, FROM the bottleneck
vector, per-sample coefficients \\\theta_k(x)\\ (length \\B\\), and the
logit is \$\$s(x) = \sum_k \theta_k(x)\\ b_k(x) + c,\$\$ with a learned
scalar bias \\c\\. The faithful self-explaining form \\\theta(x)\cdot
b(x)\\ is implemented here (the paper's contribution; a global linear
head \\s(x)=w\cdot b(x)+c\\ with constant coefficients is the simpler
defensible variant – see SEVEN_evidence.md for the flagged design
choice). The \\\theta\\-net uses `relu` / `tanh` hidden activation and a
\\\mathrm{Linear}(\cdot \to B)\\ output, with NO BatchNorm (every
operation is a function of \\x\\ alone, so the forward is per-row and
single-sample-safe).

**Torch is confined to FIT.** The bottleneck + \\\theta\\-net + bias are
trained end-to-end with torch (Adam, `BCEWithLogitsLoss`) on the rCLR'd,
universe-aligned training matrix. After training, the centred (zero-sum)
bottleneck \\W\\, the \\\theta\\-net weight matrices/biases, and the
scalar bias are EXPORTED to R as plain numeric matrices/vectors; the
python module is DISCARDED. The fitted model holds NO external pointer
and NO python dependency at score time: the log-contrast forward, the
\\\theta\\-net forward, and the logit are all pure base-R. This makes
the default `model_digest` a stable deterministic snapshot and lets the
package suite run the scorer on any host without the venv.

**Fit/score consistency despite float32.** Torch trains in float32; an R
float64 forward on the exported (float32-valued) weights differs from
the torch forward by \\\sim 10^{-6}\\. There are no frozen training
statistics beyond the exported weights themselves (the logit is a
deterministic function of the exported parameters), so fit and score
share one code path by construction – the scorer evaluates exactly the
exported parameters with `.coda_deepcoda_logit`, which is also what fit
uses to verify the held-out forward (the deepmaha lesson: fit and score
must share one code path; compute any frozen quantity from the pure-R
forward the scorer uses, NOT from torch float32 embeddings).

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
ORIGINAL abundances, NOT on `all(z == 0)`. (For a flat composition every
\\b_k = W_k \cdot 0 = 0\\, so the logit reduces to the learned bias
\\c\\, a genuine computed value, not a floored 0.)

## References

Quinn TP, Nguyen D, Rana S, Gupta S, Venkatesh S. (2020) DeepCoDA:
personalized interpretability for compositional health data.
*Proceedings of the 37th International Conference on Machine Learning
(ICML)*, PMLR 119. arXiv:2006.01392.
