# OmicSelector 2.3.0.9000 (2026-05-04, in development; v2.4.0 release accompanies Paper 3 publishable manuscript)

- Added Paper 3 within-sample compositional methods (Module A, P1) in
  `R/paper3-within-sample.R`: `ws_rclr_trimmed()` (robust trimmed CLR
  with hemolysis-marker exclusion), `ws_balance_ilr()` (orthonormal ILR
  balances on a frozen sequential binary partition encoding circulating
  miRNA biology), `ws_alr_pivot()` (additive log-ratio with a frozen
  6-miRNA pivot pool, fail-closed on missing pivots), `ws_mad_logratio()`
  (median-centred log-ratio with optional MAD scaling), and
  `ws_dominance_score()` / `fit_dominance_threshold()` /
  `ws_dominance_flag()` (panel-internal QC dominance score for
  pre-analytical-failure detection). Helpers `ws_default_sbp()` and
  `ws_default_pivot_pool()` return the curated v1 SBP tree and v1 pivot
  pool used by the balance and ALR methods respectively.
- All five methods are single-sample and reference-cohort-free: each
  sample is normalised, scored, or projected using only its own panel
  values, with no requirement for an external training cohort to exist
  at deployment time. This is the property that makes the methods
  suitable for clinical translation of circulating-miRNA biomarkers
  (Paper 3 framing: per-sample denoising for batch-effect-free
  diagnostic tests).
- Added synthetic-data unit tests in
  `tests/testthat/test-paper3-within-sample.R`.
- Pending for v2.4.0 release with Paper 3 v0.4 publishable manuscript:
  matched-null benchmark (`os_panel_null_benchmark()` extension),
  block-aware Benjamini-Hochberg with provenance blocks, MIMAT ↔
  miR-name resolver for Toray and other MIMAT-keyed deposits,
  hemolysis-injection benchmark, frozen-RUV factor estimator,
  conformal anomaly p-value with Mondrian group-conditional coverage,
  outcome-dictionary tooling, and per-cohort missingness audit. Module E
  (TabPFN, biolord, GroupDRO, horseshoe) remains scaffolded only.

# OmicSelector 2.3.0

- Added bias-audit module (6 exported functions) for cohort-provenance and
  method-ceiling diagnostics.
- Added grouped resampling helpers, provenance-floor diagnostics,
  identifiability gates, matched random-panel null benchmarks, operating-point
  summaries, and hemolysis-aware panel robustness utilities for public-omics
  biomarker-panel validation.
- Cleaned package release metadata and public documentation for the 2.3.0
  source release.

# OmicSelector 2.2.0

- Added within-sample biomarker panel methods: CLR transforms,
  CLR + MLP training and prediction, pairwise log-ratio CNN multi-seed
  reproducibility, and a CoDaCoRe interface with deterministic fallback.
- Registered new `mlr3` learners on package load: `classif.ratio_cnn`,
  `classif.clr_mlp`, and `classif.codacore`.
- Added `ws_perturbation_benchmark()` for baseline, shift, scale, and combined
  perturbation benchmarking with per-condition AUC summaries.
- Expanded test coverage for image-encoding identities, additive-shift
  invariance, CLR invariance, and learner smoke tests.
- Added CoDA-aware feature selection methods:
  `codaFS_plr_variance()`, `codaFS_selbal_wrapper()`,
  `codaFS_codacore_wrapper()`, `codaFS_logcontrast_lasso()`, and
  `codaFS_stability_logratio()`.
- Registered the new CoDA selectors as `mlr3filters` filters and wired them
  into `OmicPipeline` under the filter ids `coda_plr_variance`,
  `coda_selbal`, `coda_codacore`, `coda_logcontrast_lasso`, and
  `coda_stability_logratio`.
- Added unit tests covering additive-shift invariance, zero-sum
  log-contrast projection, synthetic discriminative log-ratio recovery, and
  smoke-scale execution at 14, 100, and 500 features.
