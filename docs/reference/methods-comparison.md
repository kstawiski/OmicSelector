# Methods Comparison Utilities

Helpers for comparing biomarker preprocessing pipelines on the same
cohort under leakage-free cross-validation. Includes:

- [`os_log_transform_adaptive()`](https://kstawiski.github.io/OmicSelector/reference/os_log_transform_adaptive.md)
  — per-feature pseudocount log2 transform that survives mixtures of raw
  counts and pre-normalized ratios.

- `os_plate_median_frozen` — R6 class for fold-frozen plate-median (or
  batch-median) centering, simpler and more deployable than ComBat when
  only a location shift is suspected.

- [`os_paired_delong()`](https://kstawiski.github.io/OmicSelector/reference/os_paired_delong.md)
  — paired DeLong test wrapper for AUC comparison of two scoring
  functions on the same cohort.

- [`os_oof_pipeline_compare()`](https://kstawiski.github.io/OmicSelector/reference/os_oof_pipeline_compare.md)
  — grouped-CV driver that runs each pipeline train-fold-frozen, pools
  out-of-fold predictions per seed, and reports paired DeLong tests vs a
  reference pipeline.

## Details

These functions exist to make the "is method X better than the current
pipeline?" question answerable in a single call instead of an ad-hoc
script. They are meant to compose with existing OmicSelector primitives
such as
[`clr_transform`](https://kstawiski.github.io/OmicSelector/reference/clr_transform.md),
[`FrozenComBat`](https://kstawiski.github.io/OmicSelector/reference/FrozenComBat.md),
[`ws_rank`](https://kstawiski.github.io/OmicSelector/reference/ws_rank.md),
[`os_panel_null_benchmark`](https://kstawiski.github.io/OmicSelector/reference/os_panel_null_benchmark.md),
and
[`os_clustered_bootstrap_auc`](https://kstawiski.github.io/OmicSelector/reference/os_clustered_bootstrap_auc.md).

## References

DeLong ER, DeLong DM, Clarke-Pearson DL. (1988). Comparing the areas
under two or more correlated receiver operating characteristic curves: a
nonparametric approach. *Biometrics*, 44(3), 837-845.

Robin X, Turck N, Hainard A, et al. (2011). pROC: an open-source package
for R and S+ to analyze and compare ROC curves. *BMC Bioinformatics*,
12, 77.

Stouffer SA, Suchman EA, DeVinney LC, Star SA, Williams RM. (1949). *The
American Soldier: Adjustment During Army Life*. Princeton University
Press.
