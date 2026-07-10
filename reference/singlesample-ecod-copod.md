# ECOD + COPOD novelty discriminator (single-sample, frozen-ECDF)

Single-sample-deployable novelty/outlier-detection-as-discriminator
scorer (family NC) that bundles the two
empirical-cumulative-distribution outlier detectors ECOD (Li et al.,
2022) and COPOD (Li et al., 2020) from the Python pyod library. Both
detectors flag a specimen as anomalous when its per-feature values fall
in the extreme tails of an in-distribution reference, which here is the
training CONTROL cohort (mirroring the healthy-reference convention of
[`fit_conformal_anomaly`](https://kstawiski.github.io/OmicSelector/reference/fit_conformal_anomaly.md)).

**Why the per-feature ECDF state is frozen rather than re-fitted at
score time.** pyod's `ECOD$decision_function()` /
`COPOD$decision_function()` are *transductive*: they re-fit each
per-feature empirical CDF on whatever batch is passed, so the score of
one specimen depends on which other specimens are scored alongside it
(measured batch-vs-row-by-row max difference \\\approx 0.68\\ for ECOD).
That is a batch-coupled scorer and would fail the single-sample
deployability gate
([`singlesample_assert_row_equivariant`](https://kstawiski.github.io/OmicSelector/reference/singlesample_assert_row_equivariant.md)).
This implementation therefore NEVER calls pyod at score time. At fit it
freezes the training-control per-feature order statistics and the
per-feature skewness signs that ECOD uses; at score it evaluates each
specimen's ECOD and COPOD aggregates against those frozen reference
ECDFs in pure R, one specimen at a time. On the training-control rows
the frozen reference ECDF equals the fit-time ECDF, so the pure-R
aggregates reproduce pyod's `decision_scores_` to machine precision
(asserted by the test suite).

**Single-sample transform.** Each specimen is mapped to the
self-contained per-sample robust CLR over its OWN strictly-positive
support (geometric-mean centring on `v > 0`); this is exactly invariant
to per-specimen scaling and uses no cross-row statistic. Missing
universe features are absent from a specimen's rCLR support and
contribute no tail term (documented single-sample-safe rule); the
per-feature reference ECDF is simply skipped for an absent feature,
exactly as if that coordinate carried the neutral aggregate `0`.

**Combining the two detectors and orientation.** ECOD and COPOD live on
different (feature-count-dependent) scales, so each per-specimen
aggregate is z-standardized against the frozen control-aggregate mean /
sd captured at fit, and the two z-scores are averaged. A single frozen
orientation scalar (\\+1\\ or \\-1\\), chosen at fit so that larger =
more case-like (training AUC of the combined anomaly score against
`y_train`; flipped if anomaly is anti-correlated with case), is applied
at score. Orientation and standardization are frozen constants, so the
score of a specimen depends only on that specimen and the frozen model.

Degenerate inputs (no present features, an all-zero specimen, or an
empty control reference) return the neutral score `0`.

## References

Li Z, Zhao Y, Botta N, Ionescu C, Hu X. (2020) COPOD: Copula-Based
Outlier Detection. *ICDM*.

Li Z, Zhao Y, Hu X, Botta N, Ionescu C, Chen G. (2022) ECOD:
Unsupervised Outlier Detection Using Empirical Cumulative Distribution
Functions. *IEEE TKDE*.

Vovk V, Gammerman A, Shafer G. (2005) Algorithmic Learning in a Random
World. Springer.
