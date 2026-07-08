# Paper 3 matched-null benchmark for panel-vs-random-panel AUC inference

Generalised matched-null benchmark introduced in Paper 3 (Module A
validation framework; Stawiski et al., in preparation). Unlike the
simpler
[`os_panel_null_benchmark`](https://kstawiski.github.io/OmicSelector/reference/os_panel_null_benchmark.md)
in `panel-gates.R`, this implementation stratifies random-panel draws by
per-feature detection-rate and log-mean-abundance quartile bins,
explicitly excludes hemolysis-marker miRNAs when
`post_hemolysis_corrected = TRUE`, and provides a three-tier fallback
protocol for small platforms with exhausted candidate pools. BH-FDR and
block-aware BH-FDR correction helpers are also provided for
multi-modality family testing.

Methods provided:

- [`paper3_matched_null_benchmark`](https://kstawiski.github.io/OmicSelector/reference/paper3_matched_null_benchmark.md):
  main benchmark function.

- [`paper3_hanley_mcneil_auc_ci`](https://kstawiski.github.io/OmicSelector/reference/paper3_hanley_mcneil_auc_ci.md):
  closed-form Hanley-McNeil 1982 95% AUC confidence interval helper.

- [`paper3_bh_fdr_correct_matched_null`](https://kstawiski.github.io/OmicSelector/reference/paper3_bh_fdr_correct_matched_null.md):
  BH-FDR within a modality family.

- [`paper3_holm_correct_familywise`](https://kstawiski.github.io/OmicSelector/reference/paper3_holm_correct_familywise.md):
  Holm correction for family-wise method × claim-state contrasts.

- [`paper3_bh_fdr_correct_blocked`](https://kstawiski.github.io/OmicSelector/reference/paper3_bh_fdr_correct_blocked.md):
  block-aware two-stage BH-FDR for specimen-shared cohort clusters.

## References

Benjamini Y, Hochberg Y. (1995) Controlling the False Discovery Rate: A
Practical and Powerful Approach to Multiple Testing. *Journal of the
Royal Statistical Society Series B* 57(1): 289–300.

Stawiski K. (in preparation) Provenance-aware within-sample scoring for
circulating-microRNA biomarkers across cancers and platforms (Paper 3 of
the OmicSelector programme; Nature Methods target).
