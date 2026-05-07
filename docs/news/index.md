# Changelog

## OmicSelector 2.4.0 (2026-05-07)

Round 3 additions (2026-05-07) — manuscript-alignment patches (Paper 3
audit patches 1 through 9; orchestrator patch 10 deferred):

- Added qPCR non-detect handling in `R/paper3-nondetects.R`:
  [`qpcr_nondetect_impute()`](https://kstawiski.github.io/OmicSelector/reference/qpcr_nondetect_impute.md)
  (Bayesian hierarchical imputation via the `nondetects` Bioconductor
  package, with fall-through to LOD when `nondetects`, `HTqPCR`, or
  `Biobase` is unavailable) and
  [`qpcr_nondetect_lod_fallback()`](https://kstawiski.github.io/OmicSelector/reference/qpcr_nondetect_lod_fallback.md)
  (limit-of-detection fallback at Ct = 40). Mirrors manuscript Methods
  §“Non-detect handling for quantitative PCR data”. Added `nondetects`,
  `HTqPCR`, and `Biobase` to `Suggests:`.
- Added
  [`hemolysis_index_blondal()`](https://kstawiski.github.io/OmicSelector/reference/hemolysis_index_blondal.md)
  in `R/paper3-hemolysis.R`: Blondal-style log(miR-451a) −
  log(miR-23a-3p) helper accepting both canonical mature miRNA names and
  MIMAT accessions. Closes the manuscript-claimed Blondal-canonical
  hemolysis index gap; pair with
  [`fit_hemolysis_rr()`](https://kstawiski.github.io/OmicSelector/reference/fit_hemolysis_rr.md)
  /
  [`apply_hemolysis_rr()`](https://kstawiski.github.io/OmicSelector/reference/apply_hemolysis_rr.md)
  for the manuscript-specified Module B nuisance correction.
- Added
  [`os_provenance_preflight()`](https://kstawiski.github.io/OmicSelector/reference/os_provenance_preflight.md)
  in `R/paper3-provenance-preflight.R`: specimen-overlap pre-flight gate
  that reads a TSV manifest and returns one of `NO_OVERLAP`,
  `KNOWN_OVERLAP`, or `UNKNOWN_ACCESSION`. Ships a curated default
  manifest at `inst/extdata/provenance_manifest.tsv` covering the
  Toray-3D cluster, Toray-V20 cluster, COSMOS-LDCT pair, VUMC-PDAC
  NanoString pair, and the audited-independent accessions used by the
  manuscript.
- Added
  [`paper3_matched_null_benchmark_cv()`](https://kstawiski.github.io/OmicSelector/reference/paper3_matched_null_benchmark_cv.md)
  in `R/paper3-matched-null.R`: the 5-fold nested-CV matched-null
  benchmark (panel selection by univariate AUC on the training fold;
  matched-null strata computed on the training fold; held-out test fold
  scored). Supports grouped folds via `group_id`, deterministic per-draw
  seeding via the `v0.6_per_draw_seed_schedule_independent` RNG
  protocol, and the three-tier matched-null fallback. This is the engine
  that produced the manuscript’s per-cell numbers; promoting it from
  `OmicSelector_paper/code/methods/matched_null_benchmark.R` makes the
  headline tables and figures reproducible from the package alone.
- Extended
  [`paper3_bh_fdr_correct_blocked()`](https://kstawiski.github.io/OmicSelector/reference/paper3_bh_fdr_correct_blocked.md):
  now reports BOTH the conservative reference (`q_block_BH`, `p_block`;
  S/c against chi-square(2k)) and the textbook Brown 1975 /
  Kost-McDermott 2002 reference (`q_block_BH_textbook`,
  `p_block_textbook`; S/c against chi-square(2k/c)). Added a `rho`
  argument (default 0.25) so the manuscript-required sensitivity
  analysis at rho ∈ {0.10, 0.50, 1.00} can be reproduced from the
  package. **Backwards-compatible:** existing call signatures continue
  to work and the `rho` default reproduces the previous numerical
  behaviour exactly. New columns are *added*, not replaced.
- Added `auc_ci_method` argument to both
  [`paper3_matched_null_benchmark()`](https://kstawiski.github.io/OmicSelector/reference/paper3_matched_null_benchmark.md)
  and
  [`paper3_matched_null_benchmark_cv()`](https://kstawiski.github.io/OmicSelector/reference/paper3_matched_null_benchmark_cv.md):
  `"hanley_mcneil"` (default; closed-form Hanley-McNeil 1982 SE),
  `"delong"` (via `pROC::ci.auc(method = "delong")`; pooled fold-level
  scores in the CV variant), or `"none"`. Returns `auc_obs_ci_lo` /
  `auc_obs_ci_hi` alongside `auc_obs` / `auc_obs_cv`. Closes the v5
  per-cell-AUC-CI gap flagged in the manuscript Discussion.
- Exported
  [`paper3_hanley_mcneil_auc_ci()`](https://kstawiski.github.io/OmicSelector/reference/paper3_hanley_mcneil_auc_ci.md)
  as the user-facing Hanley-McNeil 1982 AUC confidence-interval helper
  used internally by the matched-null benchmark functions.
- Added `min_balance_coverage` argument to
  [`ws_balance_ilr()`](https://kstawiski.github.io/OmicSelector/reference/ws_balance_ilr.md)
  (default 0.8, per manuscript Methods §“Within-sample compositional
  methods”). When the input panel computes fewer than 80% of the eight
  balances, the returned vector / matrix is filled with NA and carries
  `attr(., "coverage_failed") <- TRUE` and `attr(., "coverage")` to let
  the caller distinguish coverage failure from a real all-NA cell.
  **Behaviour change:** previously the function silently returned a
  partially-NA vector. Existing callers that need the old behaviour can
  pass `min_balance_coverage = 0`.
- Documented the
  [`ws_alr_pivot()`](https://kstawiski.github.io/OmicSelector/reference/ws_alr_pivot.md)
  fail-closed default in roxygen `@details`. Default remains
  `allow_global_fallback = FALSE`; the manuscript pipeline opts in via
  `allow_global_fallback = TRUE` with per-cell logging. No code change
  beyond the docstring.
- Bumped DESCRIPTION version from 2.3.0.9000 to 2.4.0.
- Release-check notes: `R CMD check --as-cran` may report CRAN incoming
  NOTEs for a new submission, tarball size, and optional `catboost` /
  `selbal` backends that are not in the mainstream CRAN/Bioconductor
  repositories. These are optional backend/release-distribution notes,
  not package-code warnings.
- Added unit tests for every new entry point: `tests/testthat/`
  `test-paper3-nondetects.R`, `test-paper3-provenance-preflight.R`,
  `test-paper3-hemolysis.R` (Blondal index block),
  `test-paper3-matched-null.R` (rho / textbook / per-cell AUC CI /
  nested-CV blocks), and `test-paper3-within-sample.R` (ILR coverage
  block).

Patch 10 (end-to-end pipeline orchestrator wrapping
`OmicSelector_paper/code/run_paper3_pipeline.R` as
`omicselector_paper3_run()`) is deferred to a follow-up release; it is
the only manuscript-alignment gap that remains after this round.

## OmicSelector 2.3.0.9000 (2026-05-04 / 2026-05-05, in development; v2.4.0 release accompanies Paper 3 publishable manuscript)

Round 2 additions (2026-05-05):

- Added Paper 3 miRNA name resolver in `R/paper3-mirna-name-resolver.R`:
  [`mirna_alias_table()`](https://kstawiski.github.io/OmicSelector/reference/mirna_alias_table.md)
  (88 curated miRNAs from miRBase v22.1 covering the within-sample
  sequential-binary-partition dictionary, the additive log-ratio pivot
  pool, the haemolysis markers, and Mitchell 2008 / miRBiT /
  Toray-FirePlex frequent targets),
  [`resolve_mirna_aliases()`](https://kstawiski.github.io/OmicSelector/reference/resolve_mirna_aliases.md),
  and
  [`apply_mirna_aliases()`](https://kstawiski.github.io/OmicSelector/reference/apply_mirna_aliases.md).
  Unblocks within-sample isometric-log-ratio scoring on platforms that
  deposit MIMAT or platform-specific probe IDs (Toray 3D-Gene,
  Affymetrix miRNA-3_0/-4_0, Agilent miRNA arrays) by resolving feature
  namespaces back to canonical mature-miRNA names before the
  biology-keyed partition is applied.
- Added Paper 3 matched-null benchmark in `R/paper3-matched-null.R`:
  [`paper3_matched_null_benchmark()`](https://kstawiski.github.io/OmicSelector/reference/paper3_matched_null_benchmark.md)
  (the joint detection-rate × abundance stratified matched-null
  benchmark; coexists with the simpler
  [`os_panel_null_benchmark()`](https://kstawiski.github.io/OmicSelector/reference/os_panel_null_benchmark.md)
  stub in `R/panel-gates.R`),
  [`paper3_bh_fdr_correct_matched_null()`](https://kstawiski.github.io/OmicSelector/reference/paper3_bh_fdr_correct_matched_null.md),
  [`paper3_holm_correct_familywise()`](https://kstawiski.github.io/OmicSelector/reference/paper3_holm_correct_familywise.md),
  [`paper3_bh_fdr_correct_blocked()`](https://kstawiski.github.io/OmicSelector/reference/paper3_bh_fdr_correct_blocked.md).
- Added Paper 3 preprocessing in `R/paper3-preprocessing.R`:
  [`preprocess_inverse_log()`](https://kstawiski.github.io/OmicSelector/reference/preprocess_inverse_log.md)
  (inverse-log preprocessing for pre-log-transformed microarray
  deposits).
- Added Paper 3 batch correction in `R/paper3-batch-correction.R`:
  [`fit_frozen_ruv()`](https://kstawiski.github.io/OmicSelector/reference/fit_frozen_ruv.md)
  /
  [`apply_frozen_ruv()`](https://kstawiski.github.io/OmicSelector/reference/apply_frozen_ruv.md)
  (frozen RUV factor estimator),
  [`fit_robust_pca_residual()`](https://kstawiski.github.io/OmicSelector/reference/fit_robust_pca_residual.md)
  /
  [`apply_robust_pca_residual()`](https://kstawiski.github.io/OmicSelector/reference/apply_robust_pca_residual.md).
- Added Paper 3 robust-regression haemolysis correction in
  `R/paper3-hemolysis.R`:
  [`fit_hemolysis_rr()`](https://kstawiski.github.io/OmicSelector/reference/fit_hemolysis_rr.md)
  /
  [`apply_hemolysis_rr()`](https://kstawiski.github.io/OmicSelector/reference/apply_hemolysis_rr.md)
  (per-feature M-estimator regression of expression on haemolysis index;
  coexists with the existing
  [`fit_hemolysis_prefilter()`](https://kstawiski.github.io/OmicSelector/reference/fit_hemolysis_prefilter.md)
  marker-ratio gating method in `R/hemolysis-correction.R`).
- Added Paper 3 outlier detection in `R/paper3-outlier-detection.R`:
  [`fit_compositional_mahalanobis()`](https://kstawiski.github.io/OmicSelector/reference/fit_compositional_mahalanobis.md)
  /
  [`apply_compositional_mahalanobis()`](https://kstawiski.github.io/OmicSelector/reference/apply_compositional_mahalanobis.md)
  (minimum-covariance-determinant Mahalanobis distance on log-ratio
  coordinates),
  [`fit_conformal_anomaly()`](https://kstawiski.github.io/OmicSelector/reference/fit_conformal_anomaly.md)
  /
  [`os_conformal_anomaly()`](https://kstawiski.github.io/OmicSelector/reference/os_conformal_anomaly.md)
  (conformal anomaly p-value with held-out calibration partition),
  [`fit_isolation_forest_logratio()`](https://kstawiski.github.io/OmicSelector/reference/fit_isolation_forest_logratio.md)
  /
  [`apply_isolation_forest_logratio()`](https://kstawiski.github.io/OmicSelector/reference/apply_isolation_forest_logratio.md)
  (pure-R isolation forest on rCLR inputs).
- Added Paper 3 additional within-sample methods in
  `R/paper3-additional-within-sample.R`:
  [`fit_logistic_normal_eb()`](https://kstawiski.github.io/OmicSelector/reference/fit_logistic_normal_eb.md)
  /
  [`apply_logistic_normal_eb()`](https://kstawiski.github.io/OmicSelector/reference/apply_logistic_normal_eb.md)
  (frozen-reference Empirical-Bayes denoiser),
  [`fit_frozen_quantile()`](https://kstawiski.github.io/OmicSelector/reference/fit_frozen_quantile.md)
  /
  [`apply_frozen_quantile()`](https://kstawiski.github.io/OmicSelector/reference/apply_frozen_quantile.md)
  (monotone quantile calibrator for cross-platform mapping).
- Added 6 new test files in `tests/testthat/` covering 44 test_that
  blocks for the round-2 methods (1 skip on
  `fit_compositional_mahalanobis` headline test when `robustbase` not
  installed, handled by `skip_if_not_installed`). All other tests pass.

Round 1 additions (2026-05-04):

- Added Paper 3 within-sample compositional methods (Module A, P1) in
  `R/paper3-within-sample.R`:
  [`ws_rclr_trimmed()`](https://kstawiski.github.io/OmicSelector/reference/ws_rclr_trimmed.md)
  (robust trimmed CLR with hemolysis-marker exclusion),
  [`ws_balance_ilr()`](https://kstawiski.github.io/OmicSelector/reference/ws_balance_ilr.md)
  (orthonormal ILR balances on a frozen sequential binary partition
  encoding circulating miRNA biology),
  [`ws_alr_pivot()`](https://kstawiski.github.io/OmicSelector/reference/ws_alr_pivot.md)
  (additive log-ratio with a frozen 6-miRNA pivot pool, fail-closed on
  missing pivots),
  [`ws_mad_logratio()`](https://kstawiski.github.io/OmicSelector/reference/ws_mad_logratio.md)
  (median-centred log-ratio with optional MAD scaling), and
  [`ws_dominance_score()`](https://kstawiski.github.io/OmicSelector/reference/ws_dominance_score.md)
  /
  [`fit_dominance_threshold()`](https://kstawiski.github.io/OmicSelector/reference/fit_dominance_threshold.md)
  /
  [`ws_dominance_flag()`](https://kstawiski.github.io/OmicSelector/reference/ws_dominance_flag.md)
  (panel-internal QC dominance score for pre-analytical-failure
  detection). Helpers
  [`ws_default_sbp()`](https://kstawiski.github.io/OmicSelector/reference/ws_default_sbp.md)
  and
  [`ws_default_pivot_pool()`](https://kstawiski.github.io/OmicSelector/reference/ws_default_pivot_pool.md)
  return the curated v1 SBP tree and v1 pivot pool used by the balance
  and ALR methods respectively.
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
  matched-null benchmark
  ([`os_panel_null_benchmark()`](https://kstawiski.github.io/OmicSelector/reference/os_panel_null_benchmark.md)
  extension), block-aware Benjamini-Hochberg with provenance blocks,
  MIMAT ↔︎ miR-name resolver for Toray and other MIMAT-keyed deposits,
  hemolysis-injection benchmark, frozen-RUV factor estimator, conformal
  anomaly p-value with Mondrian group-conditional coverage,
  outcome-dictionary tooling, and per-cohort missingness audit. Module E
  (TabPFN, biolord, GroupDRO, horseshoe) remains scaffolded only.

## OmicSelector 2.3.0

- Added bias-audit module (6 exported functions) for cohort-provenance
  and method-ceiling diagnostics.
- Added grouped resampling helpers, provenance-floor diagnostics,
  identifiability gates, matched random-panel null benchmarks,
  operating-point summaries, and hemolysis-aware panel robustness
  utilities for public-omics biomarker-panel validation.
- Cleaned package release metadata and public documentation for the
  2.3.0 source release.

## OmicSelector 2.2.0

- Added within-sample biomarker panel methods: CLR transforms, CLR + MLP
  training and prediction, pairwise log-ratio CNN multi-seed
  reproducibility, and a CoDaCoRe interface with deterministic fallback.
- Registered new `mlr3` learners on package load: `classif.ratio_cnn`,
  `classif.clr_mlp`, and `classif.codacore`.
- Added
  [`ws_perturbation_benchmark()`](https://kstawiski.github.io/OmicSelector/reference/ws_perturbation_benchmark.md)
  for baseline, shift, scale, and combined perturbation benchmarking
  with per-condition AUC summaries.
- Expanded test coverage for image-encoding identities, additive-shift
  invariance, CLR invariance, and learner smoke tests.
- Added CoDA-aware feature selection methods:
  [`codaFS_plr_variance()`](https://kstawiski.github.io/OmicSelector/reference/codaFS_plr_variance.md),
  [`codaFS_selbal_wrapper()`](https://kstawiski.github.io/OmicSelector/reference/codaFS_selbal_wrapper.md),
  [`codaFS_codacore_wrapper()`](https://kstawiski.github.io/OmicSelector/reference/codaFS_codacore_wrapper.md),
  [`codaFS_logcontrast_lasso()`](https://kstawiski.github.io/OmicSelector/reference/codaFS_logcontrast_lasso.md),
  and
  [`codaFS_stability_logratio()`](https://kstawiski.github.io/OmicSelector/reference/codaFS_stability_logratio.md).
- Registered the new CoDA selectors as `mlr3filters` filters and wired
  them into `OmicPipeline` under the filter ids `coda_plr_variance`,
  `coda_selbal`, `coda_codacore`, `coda_logcontrast_lasso`, and
  `coda_stability_logratio`.
- Added unit tests covering additive-shift invariance, zero-sum
  log-contrast projection, synthetic discriminative log-ratio recovery,
  and smoke-scale execution at 14, 100, and 500 features.
