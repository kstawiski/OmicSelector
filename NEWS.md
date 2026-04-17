# OmicSelector 2.2.0

- Added the Paper 1 v2.2 within-sample biomarker panel methods: CLR transforms,
  CLR + MLP training and prediction, pairwise log-ratio CNN multi-seed
  reproducibility, and a CoDaCoRe interface with deterministic fallback.
- Registered new `mlr3` learners on package load: `classif.ratio_cnn`,
  `classif.clr_mlp`, and `classif.codacore`.
- Added `ws_perturbation_benchmark()` for baseline, shift, scale, and combined
  perturbation benchmarking with per-condition AUC summaries.
- Expanded test coverage for image-encoding identities, additive-shift
  invariance, CLR invariance, and learner smoke tests.
- Added CoDA-aware feature selection methods for Paper 1 v2.3:
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
