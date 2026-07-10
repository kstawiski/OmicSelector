# Single-sample matched-null benchmark for panel-vs-random-panel AUC inference

Generalised matched-null benchmark introduced in the single-sample
scoring bank (Module A validation framework). Unlike the simpler
[`os_panel_null_benchmark`](https://kstawiski.github.io/OmicSelector/reference/os_panel_null_benchmark.md)
in `panel-gates.R`, this implementation stratifies random-panel draws by
per-feature detection-rate and log-mean-abundance quartile bins,
explicitly excludes hemolysis-marker miRNAs when
`post_hemolysis_corrected = TRUE`, and provides a three-tier fallback
protocol for small platforms with exhausted candidate pools. BH-FDR and
block-aware BH-FDR correction helpers are also provided for
multi-modality family testing.

Methods provided:

- [`singlesample_matched_null_benchmark`](https://kstawiski.github.io/OmicSelector/reference/singlesample_matched_null_benchmark.md):
  main benchmark function.

- [`singlesample_hanley_mcneil_auc_ci`](https://kstawiski.github.io/OmicSelector/reference/singlesample_hanley_mcneil_auc_ci.md):
  closed-form Hanley-McNeil 1982 95% AUC confidence interval helper.

- [`singlesample_bh_fdr_correct_matched_null`](https://kstawiski.github.io/OmicSelector/reference/singlesample_bh_fdr_correct_matched_null.md):
  BH-FDR within a modality family.

- [`singlesample_holm_correct_familywise`](https://kstawiski.github.io/OmicSelector/reference/singlesample_holm_correct_familywise.md):
  Holm correction for family-wise method × claim-state contrasts.

- [`singlesample_bh_fdr_correct_blocked`](https://kstawiski.github.io/OmicSelector/reference/singlesample_bh_fdr_correct_blocked.md):
  block-aware two-stage BH-FDR for specimen-shared cohort clusters.

## References

Benjamini Y, Hochberg Y. (1995) Controlling the False Discovery Rate: A
Practical and Powerful Approach to Multiple Testing. *Journal of the
Royal Statistical Society Series B* 57(1): 289–300.
