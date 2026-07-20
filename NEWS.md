# OmicSelector 2.6.5.9001 (development)

- Cache each specimen's detrended MFDFA window variances once per scale in
  `frac-mfdfa` and reuse them across q values. This removes sixfold redundant
  least-squares work under the default q grid while preserving the legacy
  descriptor exactly.

- Clarify the paper benchmark's held-out support contract so complete finite
  constant scores are retained as valid midrank-AUC evidence (AUC 0.5), while
  incomplete or non-finite score vectors still fail closed. The external
  rank-tree runner now applies this rule uniformly, binds both approved design
  amendments, and resolves its installed dependency closure only within the
  active R libraries so a pinned `pkgload` source checkout cannot invalidate
  final runtime re-attestation. The shared paper runtime-receipt producer now
  uses the same installed-library resolution contract.

- Correct `os_detect_cross_cohort_duplicates()` so its documented `tolerance`
  and `min_features` arguments are enforced. Exact matching no longer rounds
  profiles to eight decimals, tolerant/partially observed matching is explicit,
  and malformed identifiers or thresholds fail closed.
- Bound the external `ranktreeEnsemble` adapters' training-only preliminary
  importance screen to a fixed 128-feature target before quadratic pair
  expansion. This preserves the backend's top-quartile default on ordinary
  panels while preventing protection/C-stack failure on high-dimensional miRNA
  inputs. Prediction now reconstructs and requires only the exact screened
  raw-feature subset rather than requiring and expanding all original inputs
  again; the frozen target, realized backend quantile, and selected feature
  order are stored in every model.
- Add audited `ranktreeEnsemble` 0.24 adapters for random rank forests and
  boosted rank trees, plus an explicit external-competitor registry. These
  post-freeze sensitivity routes preserve the 74-row primary roster, construct
  every learned representation inside training folds, enforce a frozen feature
  map and use `datrank = FALSE` to avoid the backend's train/predict
  representation mismatch. The adapter also confines random-forest screening
  to its reproducible single-core RNG route and preserves singleton probability
  array dimensions.
- Add audited group-collapsed and profile-weighted inference helpers for fair
  comparison of the full registered single-sample method roster, including
  corrected repeated-CV uncertainty, frozen multiplicity families, and
  outcome-independent complete-support panels. The paper runner uses the
  group-collapsed estimand as primary when biological groups have repeated
  profiles and reports equal-profile weighting as a sensitivity analysis.
- Add an experimental grouped cross-fitted single-sample selector with
  best-AUC, lower-confidence-limit, and simplex-stack routes. The three routes
  can share one inner cross-fit pass and return pure-R frozen scoring objects.
  Deployed selector scores now route every base method and the final simplex
  affine reduction through the same one-specimen arithmetic path, guaranteeing
  exact batch-versus-singleton equality across BLAS implementations.
- Add immutable paper-analysis producers for snapshot integration, all-method
  comparison, fully nested selector tasks, and selector synthesis. These
  scripts fail closed on provenance, split, support, fit-budget, and manifest
  inconsistencies.
- Add a fail-closed 34-by-five nested-selector task-grid producer and portable
  one-task launcher that pin cache bytes, clean package commit, runtime inputs,
  and atomic output paths.
- Pin nested-selector tasks to exact container and R-library snapshot hashes,
  record the executing R/package environment, and carry the frozen
  modality/biospecimen and 33-primary/one-sensitivity unit accounting into
  synthesis.
- Add a fail-closed, serialization-only splitter for the paper's frozen
  optional-dependency benchmark. It verifies the complete 13-by-five-by-34
  task contract, provenance-group folds, result cardinality, and package/code
  pins before atomically emitting one-cell task bundles and a final manifest;
  arbitrary fit, score, dependency, and runtime failures cannot be relabelled
  as structural ineligibility.
- Add a package-owned structural-witness producer for generic ineligibility in
  that dependency benchmark. It replays only explicitly pinned cells from a
  clean historical package commit, preserves every finite held-out score, and
  authorizes sharding only when exact cache, task, analysis, fold-engine, and
  producer bytes prove a fold-local constant-score limitation rather than a
  hidden fit, scoring, dependency, or runtime error.
- Make the paper's all-method comparison consume an immutable, fully validated
  dependency-overlay candidate rather than mutable aggregate tables. The
  snapshot now binds every upstream validation row and input byte, exact source
  and installed-package pins, route-specific execution receipts, splitter
  manifests, and the structural witness transitively through the candidate
  manifest. Snapshot construction, execution, and validation additionally bind
  a deterministic digest of every file in the loaded OmicSelector installation,
  preventing a same-version but stale lazy-load database from passing a
  source-script check. Fair-inference snapshots export predictions only for
  method--unit pairs completing all five registered seeds, while preserving
  every partial seed failure and reason in the complete cell/eligibility audit.
- Add a package-owned, atomic runtime-receipt producer and shared start/end
  guard for the paper analyses. It pins the immutable image, canonical
  executables, exact R/platform/RNG/BLAS and one-thread state, and recursive
  installed dependency trees rooted at OmicSelector, data.table, digest, and
  metafor; snapshot, runner, and validator manifests carry the same receipt
  and image digests.
- Add a package-owned post-run validator for the fair all-method comparison. It
  verifies the exact artifact inventory, frozen dimensions and pair family,
  common-support held-out joins, multiplicity families, inference and
  resampling contracts, provenance identity, and source/snapshot/runtime pins
  before atomically promoting a PASS receipt.

# OmicSelector 2.6.5 (2026-07-17)

* Reimplemented the exact `os_ktsp_fit()` candidate-pair enumeration as
  vectorized feature blocks with top-k retention. The fitted pair ordering and
  scores are identical to the prior pairwise implementation, while large miRNA
  feature spaces no longer incur millions of interpreted-R loop iterations.

# OmicSelector 2.6.4 (2026-07-16)

* `score_lrt_bw()` now reuses the fit-time Bures-Wasserstein class
  representation when the complete fitted feature universe is present. This
  removes amplified eigensolver round-off differences between singleton and
  batch scoring while preserving the documented partial-overlap route.

# OmicSelector 2.6.3 (2026-07-16)

- Promoted the benchmark's trimmed-rCLR signed-panel reference to canonical
  `fit_ws_rclr_panel()` / `score_ws_rclr_panel()` package functions. The frozen
  model selects a 20-feature panel and score signs from training data only and
  exposes the same one-score-per-specimen deployment contract as the other
  rostered methods.

# OmicSelector 2.6.2 (2026-07-15)

- Made `orient = "fixed"` the default for
  `singlesample_paired_auc_diff_se()`. Outcome-dependent `"median"` and
  sign-invariant `"auc"` orientation remain explicit opt-ins.

# OmicSelector 2.6.1 (2026-07-15)

- Added `orient = "fixed"` to `singlesample_paired_auc_diff_se()` so a
  training-frozen score direction can be evaluated without using held-out
  outcomes to reflect the ROC curve.

# OmicSelector 2.6.0 (2026-07-08)

## Single-sample deployment API

* Added a single-sample deployment API — `deploy_singlesample()` freezes one
  rostered scorer, `score_specimen()` scores an incoming specimen from its own
  panel values plus frozen fit-time parameters (no test batch, reference
  cohort, or batch-correction step), and `is_singlesample_deployable()` checks
  the singleton-equals-batch guarantee; about 30 roster methods deploy and the
  rest reject cleanly. This exposes deployability only and is not a
  superiority claim over batch-corrected pipelines. See
  `vignette("single-sample-deployment")`.
* inv-scatter: replaced the O(D^3) frozen ridge head with an exact SVD
  row-space solve so the scattering scorer is shippable when D >> n; results
  are unchanged up to numerical-rank tolerance (added
  `test-inv-scatter-head-equivalence.R`).

## Single-sample miRNA scoring method bank (74 methods)

* Completed the frozen 74-method single-sample scoring bank under
  `R/singlesample-*.R`, rostered in
  `inst/extdata/singlesample_method_manifest.csv` and enumerable via
  `singlesample_method_roster()` / `singlesample_method_bank()`. Every method
  is a single-sample (per-row) scorer: a deployed sample is scored using only
  its own panel values plus frozen parameters learned at fit time, with no
  batch-correction step and no requirement for a co-resident reference cohort
  at deployment. Methods carry a `within` (within-cohort) or `transfer`
  (cross-cohort generalisation) estimand and a role of baseline,
  discriminator, or negative control. The bank is grouped into the following
  families (manifest `family` codes A-O, NC, which are short grouping letters,
  not descriptive names):

* Compositional within-sample baselines (family A): trimmed robust-CLR
  reference (`ws_rclr_trimmed`, the reference, not a candidate), additive
  log-ratio with a frozen pivot pool, MAD log-ratio, and a panel-internal
  dominance score for pre-analytical-failure detection.

* Frozen denoisers (family B): frozen-RUV factor removal, robust-PCA residual,
  and robust-regression hemolysis residual removal.

* Kit-aware compositional transfer scorers (family D): kit-stratified rCLR,
  kit fixed-effect-adjusted ALR, kit-orthogonal ILR, and kit-residual MAD.

* Anchor/reference denoising transfer scorers (family E): kit-stable-anchor
  rCLR, hemolysis-and-kit dual-anchor, and RIN-weighted reference scoring.

* Technology-aware transfer scorers (family F): technology-stratified rCLR,
  cross-technology harmonised rCLR, and technology-residualised ALR.

* Biofluid-aware transfer scorers (family G): biofluid-stratified rCLR,
  biofluid-anchor rCLR, biofluid-residualised ALR, and blood-cell/platelet
  marker-exclusion rCLR.

* Rank / relative-ordering and balance methods (family K): k-TSP relative
  expression ordering (`reo_ktsp`), pairwise-log-ratio penalised logistic
  (`reo_pairratio`), singscore and UCell within-sample rank scores, a
  meta-analytic cross-cohort k-TSP (`reo_metaktsp`), the `selbal` balance
  selection + log-contrast logistic scorer (`bal_selbal`), and the ILR
  balance scorer (`ws_balance_ilr`).

* Likelihood-ratio, optimal-transport, kernel and density-ratio scorers
  (family L): class-conditional likelihood-ratio tests including
  Bures-Wasserstein (`lrt_bw`), KDE naive-Bayes (`lrt_nbkde`),
  Gaussian-copula (`lrt_copula`), tail-aware Student-t-copula (`lrt_tcopula`),
  vine-copula (`lrt_vinecopula`, via `rvinecopulib`), Mondrian-conformal LRT
  (`conf_mondrian`), linearised OT embedding to a frozen reference (`ot_lot`),
  sliced-Wasserstein-view Gaussian LRT (`ot_slicedlrt`), kernel
  mean-embedding witness (`kme_witness`), Fisher-Rao geodesic LRT
  (`ig_fisherrao`), and direct density-ratio uLSIF (`dre_ulsif`).

* Image / invariance / signal-descriptor scorers (family M): wavelet
  scattering + frozen head (`inv_scatter`, via `kymatio`), ordinal-LBP
  (`inv_olbp`), sublevel persistent homology persistence image (`tda_ph`),
  path-signature features (`sig_path`), multifractal-DFA spectrum descriptors
  (`frac_mfdfa`), GASF-image-to-frozen-CNN embedding (`img_gasfcnn`, via
  torch/torchvision), graph-Fourier transform on a frozen co-expression graph
  (`gsp_gft`), Nystrom out-of-sample diffusion-map (`man_nystrom`), Haralick
  GLCM (`inv_glcm`), quantile-function FDA (`inv_fdaqf`), and
  curvature-scale-space descriptor (`inv_css`).

* Neural / domain-generalisation / self-supervised scorers (family N): most
  use a torch backend at fit and export frozen weights for a pure-R per-row
  score. Members are domain-adversarial encoder (`dann`), counterfactual
  class-conditional VAE (`cvae`), Invariant Causal Prediction (`icp`),
  SCARF contrastive (`ai_scarf`), group-DRO (`dro_group`), V-REx
  (`dro_vrex`), self-gated mixture of frozen experts (`moe_gated`),
  spectral-normalised neural GP (`unc_sngp`), VICReg self-supervised
  embedding (`ssl_vicreg`), prototypical network (`proto_net`), deep
  Mahalanobis on a learned frozen embedding (`lrt_deepmaha`), Fishr
  (`dg_fishr`), IB-IRM (`dg_ibirm`), and StableMate stable-predictor
  selection (`sel_stablemate`).

* Foundation-model and non-linear-compositional scorers (family O):
  in-context tabular foundation models scored against a frozen reference
  context (TabPFN-v2 `tabpfn`, large-context TabICL `tabicl`, retrieval-based
  TabDPT `tabdpt`), non-linear log-contrast network DeepCoDA
  (`coda_deepcoda`), and sparse log-ratio balances CoDaCoRe (`coda_codacore`).

* Negative controls (family NC): deliberately weak atypicality / outlier
  detectors used as a calibration floor, not biomarker candidates -- frozen
  MCD Mahalanobis (`fit_compositional_mahalanobis`), conformal anomaly
  p-value (`fit_conformal_anomaly`), log-ratio isolation forest
  (`fit_isolation_forest_logratio`), ECOD/COPOD global tail outlier
  (`ecod_copod`, via PyOD), and a forced n=1 entropic-kernel single-sample
  Sinkhorn control (`sinkhorn_single`).

* Methods using a Python backend (torch, kymatio, torchvision, PyOD, TabPFN,
  TabICL, TabDPT, via `reticulate`) train at fit time and, where applicable,
  export frozen weights so the deployment-time score runs in base R; their
  optional backends are not declared as hard dependencies. The bank can be
  audited against exported package functions with
  `singlesample_assert_method_bank_exports()`.

* Vectorized the `ws_balance_ilr()` matrix/data.frame path for deployment
  throughput (a large speed-up on realistic panels) while preserving
  byte-identical results versus the prior per-row path on realistic inputs.
  This incidentally fixes data.frame recursion and single-column matrix
  name-dropping. The pre-existing `TOP_K_BY_ABUNDANCE` `p <= k` empty-tail
  quirk is intentionally preserved for frozen within-sample result identity
  and is flagged for a future separate fix because changing it alters scores.

## Maintenance

* Corrected the compositional Mahalanobis scorer to use full-rank coordinates
  rather than a singular full-composition covariance representation.
* Preserved exact row-wise TabICL scoring for large feature panels and ensured
  TabDPT receives writable, owned NumPy arrays from R inputs.
* Prevented pair-ratio logistic fits from using a truncated regularization
  path when selecting the frozen deployment model.
* Added a project `renv.lock` and activation bootstrap for reproducible package
  development and checks.
* Renamed the single-sample-deployable method-bank machinery from the internal
  codename `paper3_*` / `R/paper3-*.R` to the descriptive `singlesample_*` /
  `R/singlesample-*.R` (e.g. `singlesample_method_roster()`,
  `singlesample_score_call()`, `singlesample_assert_row_equivariant()`;
  `inst/extdata/singlesample_method_manifest.csv`). The former `paper3_*`
  exported names are retained as thin, exported back-compatibility aliases, so
  existing code keeps working; prefer the `singlesample_*` names.

# OmicSelector 2.5.0 (2026-05-29)

* Aligned promoted Paper 3 matched-null nested-CV scoring with the CLOSED A2
  F4 repair: the CV p-value now requires both per-fold null validity
  (`fold_fallback_frac < 0.5`) and at least `K/2` complete null draws across
  null-valid folds. Degenerate matched-null cells now fail closed with
  `eligible = FALSE` and `p_emp_cv = NA_real_`.
* Added `singlesample_paired_auc_diff_se()` and
  `singlesample_technology_lift_delong()` for empirical paired DeLong AUC-lift
  standard errors, matching the manuscript's F6 mean-of-folds estimand and
  technology-transfer orientation.
* Package-ized the promoted Paper 3 runner utilities: hardcoded
  manuscript-workspace input/output paths were replaced with caller-supplied
  input paths, package `extdata` lookup where available, and `tempdir()`
  output defaults.
* Verified the F5 transfer-LOPBO and F3 brms-retirement fixes are
  manuscript-driver concerns with no clean package counterpart in the exported
  Paper 3 scorer primitives; no brms/Stan alternative-inference branch is
  present in the package.

# OmicSelector 2.4.0 (2026-05-07)

Round 5 additions (2026-05-22):

* Promoted the remaining OmicSelector Paper 3 method-family scorer
  implementations from the manuscript analysis workspace into the package:
  provenance-aware, kit-aware, anchor/reference, technology-aware,
  biofluid-aware, learned kit-aware, Group-DRO, and Sinkhorn-OT scorers.
* Added `singlesample_method_bank()` and `singlesample_assert_method_bank_exports()` so the
  manuscript method bank can be audited directly against exported package
  functions.
* Added `singlesample_make_loco_splits()` and `singlesample_make_locto_splits()` for
  same-provenance-block-excluding leave-one-cohort-out and
  leave-one-cancer-type-out validation protocols.
* Added focused unit tests for the promoted method families and shipped the
  optional Sinkhorn Python helper under `inst/python/`.

Round 4 additions (2026-05-08):

* Added MIMAT ↔ hsa-miR namespace lookup table (miRBase 22.1; 2,656 human mature miRNA entries) and `resolve_mimat_to_hsa()` / `resolve_hsa_to_mimat()` helpers (`mimat_hsa_lookup_v22_1.tsv` shipped under `inst/extdata/`). Matching is case-insensitive; canonical-case output (MIMAT uppercase, hsa-miR lowercase).

Round 3 additions (2026-05-07) — manuscript-alignment patches (Paper 3 audit
patches 1 through 9; orchestrator patch 10 deferred):

- Added qPCR non-detect handling in `R/singlesample-nondetects.R`:
  `qpcr_nondetect_impute()` (Bayesian hierarchical imputation via the
  `nondetects` Bioconductor package, with fall-through to LOD when
  `nondetects`, `HTqPCR`, or `Biobase` is unavailable) and
  `qpcr_nondetect_lod_fallback()` (limit-of-detection fallback at Ct = 40).
  Mirrors manuscript Methods §"Non-detect handling for quantitative PCR
  data". Added `nondetects`, `HTqPCR`, and `Biobase` to `Suggests:`.
- Added `hemolysis_index_blondal()` in `R/singlesample-hemolysis.R`: Blondal-style
  log(miR-451a) − log(miR-23a-3p) helper accepting both canonical mature
  miRNA names and MIMAT accessions. Closes the manuscript-claimed
  Blondal-canonical hemolysis index gap; pair with `fit_hemolysis_rr()` /
  `apply_hemolysis_rr()` for the manuscript-specified Module B nuisance
  correction.
- Added `os_provenance_preflight()` in `R/singlesample-provenance-preflight.R`:
  specimen-overlap pre-flight gate that reads a TSV manifest and returns
  one of `NO_OVERLAP`, `KNOWN_OVERLAP`, or `UNKNOWN_ACCESSION`. Ships a
  curated default manifest at `inst/extdata/provenance_manifest.tsv`
  covering the Toray-3D cluster, Toray-V20 cluster, COSMOS-LDCT pair,
  VUMC-PDAC NanoString pair, and the audited-independent accessions used
  by the manuscript.
- Added `singlesample_matched_null_benchmark_cv()` in `R/singlesample-matched-null.R`:
  the 5-fold nested-CV matched-null benchmark (panel selection by
  univariate AUC on the training fold; matched-null strata computed on
  the training fold; held-out test fold scored). Supports grouped folds
  via `group_id`, deterministic per-draw seeding via the
  `v0.6_per_draw_seed_schedule_independent` RNG protocol, and the
  three-tier matched-null fallback. This is the engine that produced
  the manuscript's per-cell numbers; promoting it from
  `OmicSelector_paper/code/methods/matched_null_benchmark.R` makes the
  headline tables and figures reproducible from the package alone.
- Extended `singlesample_bh_fdr_correct_blocked()`: now reports BOTH the
  conservative reference (`q_block_BH`, `p_block`; S/c against
  chi-square(2k)) and the textbook Brown 1975 / Kost-McDermott 2002
  reference (`q_block_BH_textbook`, `p_block_textbook`; S/c against
  chi-square(2k/c)). Added a `rho` argument (default 0.25) so the
  manuscript-required sensitivity analysis at rho ∈ {0.10, 0.50, 1.00}
  can be reproduced from the package. **Backwards-compatible:** existing
  call signatures continue to work and the `rho` default reproduces
  the previous numerical behaviour exactly. New columns are *added*,
  not replaced.
- Added `auc_ci_method` argument to both `singlesample_matched_null_benchmark()`
  and `singlesample_matched_null_benchmark_cv()`: `"hanley_mcneil"` (default;
  closed-form Hanley-McNeil 1982 SE), `"delong"` (via
  `pROC::ci.auc(method = "delong")`; pooled fold-level scores in the CV
  variant), or `"none"`. Returns `auc_obs_ci_lo` / `auc_obs_ci_hi`
  alongside `auc_obs` / `auc_obs_cv`. Closes the v5 per-cell-AUC-CI gap
  flagged in the manuscript Discussion.
- Exported `singlesample_hanley_mcneil_auc_ci()` as the user-facing
  Hanley-McNeil 1982 AUC confidence-interval helper used internally by the
  matched-null benchmark functions.
- Added `min_balance_coverage` argument to `ws_balance_ilr()` (default
  0.8, per manuscript Methods §"Within-sample compositional methods").
  When the input panel computes fewer than 80% of the eight balances,
  the returned vector / matrix is filled with NA and carries
  `attr(., "coverage_failed") <- TRUE` and `attr(., "coverage")` to
  let the caller distinguish coverage failure from a real all-NA cell.
  **Behaviour change:** previously the function silently returned a
  partially-NA vector. Existing callers that need the old behaviour can
  pass `min_balance_coverage = 0`.
- Documented the `ws_alr_pivot()` fail-closed default in roxygen
  `@details`. Default remains `allow_global_fallback = FALSE`; the
  manuscript pipeline opts in via `allow_global_fallback = TRUE` with
  per-cell logging. No code change beyond the docstring.
- Bumped DESCRIPTION version from 2.3.0.9000 to 2.4.0.
- Release-check notes: `R CMD check --as-cran` may report CRAN incoming
  NOTEs for a new submission, tarball size, and optional `catboost` /
  `selbal` backends that are not in the mainstream CRAN/Bioconductor
  repositories. These are optional backend/release-distribution notes, not
  package-code warnings.
- Added unit tests for every new entry point: `tests/testthat/`
  `test-singlesample-nondetects.R`, `test-singlesample-provenance-preflight.R`,
  `test-singlesample-hemolysis.R` (Blondal index block), `test-singlesample-matched-null.R`
  (rho / textbook / per-cell AUC CI / nested-CV blocks), and
  `test-singlesample-within-sample.R` (ILR coverage block).

Patch 10 (end-to-end pipeline orchestrator wrapping
`OmicSelector_paper/code/run_singlesample_pipeline.R` as
`omicselector_singlesample_run()`) is deferred to a follow-up release; it is the
only manuscript-alignment gap that remains after this round.

# OmicSelector 2.3.0.9000 (2026-05-04 / 2026-05-05, in development; v2.4.0 release accompanies Paper 3 publishable manuscript)

Round 2 additions (2026-05-05):

- Added Paper 3 miRNA name resolver in `R/singlesample-mirna-name-resolver.R`:
  `mirna_alias_table()` (88 curated miRNAs from miRBase v22.1 covering
  the within-sample sequential-binary-partition dictionary, the additive
  log-ratio pivot pool, the haemolysis markers, and Mitchell 2008 /
  miRBiT / Toray-FirePlex frequent targets), `resolve_mirna_aliases()`,
  and `apply_mirna_aliases()`. Unblocks within-sample isometric-log-ratio
  scoring on platforms that deposit MIMAT or platform-specific probe IDs
  (Toray 3D-Gene, Affymetrix miRNA-3_0/-4_0, Agilent miRNA arrays) by
  resolving feature namespaces back to canonical mature-miRNA names
  before the biology-keyed partition is applied.
- Added Paper 3 matched-null benchmark in `R/singlesample-matched-null.R`:
  `singlesample_matched_null_benchmark()` (the joint detection-rate × abundance
  stratified matched-null benchmark; coexists with the simpler
  `os_panel_null_benchmark()` stub in `R/panel-gates.R`),
  `singlesample_bh_fdr_correct_matched_null()`,
  `singlesample_holm_correct_familywise()`,
  `singlesample_bh_fdr_correct_blocked()`.
- Added Paper 3 preprocessing in `R/singlesample-preprocessing.R`:
  `preprocess_inverse_log()` (inverse-log preprocessing for
  pre-log-transformed microarray deposits).
- Added Paper 3 batch correction in `R/singlesample-batch-correction.R`:
  `fit_frozen_ruv()` / `apply_frozen_ruv()` (frozen RUV factor estimator),
  `fit_robust_pca_residual()` / `apply_robust_pca_residual()`.
- Added Paper 3 robust-regression haemolysis correction in
  `R/singlesample-hemolysis.R`: `fit_hemolysis_rr()` / `apply_hemolysis_rr()`
  (per-feature M-estimator regression of expression on haemolysis index;
  coexists with the existing `fit_hemolysis_prefilter()` marker-ratio
  gating method in `R/hemolysis-correction.R`).
- Added Paper 3 outlier detection in `R/singlesample-outlier-detection.R`:
  `fit_compositional_mahalanobis()` / `apply_compositional_mahalanobis()`
  (minimum-covariance-determinant Mahalanobis distance on log-ratio
  coordinates), `fit_conformal_anomaly()` / `os_conformal_anomaly()`
  (conformal anomaly p-value with held-out calibration partition),
  `fit_isolation_forest_logratio()` /
  `apply_isolation_forest_logratio()` (pure-R isolation forest on rCLR
  inputs).
- Added Paper 3 additional within-sample methods in
  `R/singlesample-additional-within-sample.R`: `fit_logistic_normal_eb()` /
  `apply_logistic_normal_eb()` (frozen-reference Empirical-Bayes
  denoiser), `fit_frozen_quantile()` / `apply_frozen_quantile()`
  (monotone quantile calibrator for cross-platform mapping).
- Added 6 new test files in `tests/testthat/` covering 44 test_that
  blocks for the round-2 methods (1 skip on `fit_compositional_mahalanobis`
  headline test when `robustbase` not installed, handled by
  `skip_if_not_installed`). All other tests pass.

Round 1 additions (2026-05-04):


- Added Paper 3 within-sample compositional methods (Module A, P1) in
  `R/singlesample-within-sample.R`: `ws_rclr_trimmed()` (robust trimmed CLR
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
  `tests/testthat/test-singlesample-within-sample.R`.
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
