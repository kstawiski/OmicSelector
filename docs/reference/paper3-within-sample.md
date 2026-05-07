# Paper 3 within-sample compositional methods (Module A, P1)

Five Module-A within-sample compositional transforms introduced in Paper
3 of the OmicSelector programme (Stawiski et al., in preparation; Nature
Methods target). All methods are single-sample and
reference-cohort-free: each sample is normalised, scored, or projected
using only its own panel values, with no requirement for an external
training cohort to exist at deployment time. This is the property that
makes within-sample methods suitable for clinical translation: a
diagnostic test built on these primitives can be deployed by computing a
per-patient score from the patient's own panel measurement, without
re-fitting any normalization on a population reference.

Methods provided:

- [`ws_rclr_trimmed`](https://kstawiski.github.io/OmicSelector/reference/ws_rclr_trimmed.md):
  robust trimmed centred log-ratio with optional zero-detection-only
  convention; tolerant to dominant contaminating species (haemolysis
  miRNAs) and to small-head-of-pool compositions typical of plasma/serum
  miRNA panels.

- [`ws_balance_ilr`](https://kstawiski.github.io/OmicSelector/reference/ws_balance_ilr.md):
  orthonormal isometric log-ratio balances on a frozen sequential binary
  partition encoding circulating miRNA biology (red-blood-cell vs rest,
  platelet vs rest, let-7 family within-cluster, miR-17~92, miR-200,
  miR-371~373, dominant-vs-tail abundance balance).

- [`ws_alr_pivot`](https://kstawiski.github.io/OmicSelector/reference/ws_alr_pivot.md):
  additive log-ratio with a frozen pivot pool of putatively low-variance
  circulating miRNAs serving as the geometric-mean denominator (anchors
  the transformation against a stable internal reference rather than
  choosing a single feature as denominator).

- [`ws_mad_logratio`](https://kstawiski.github.io/OmicSelector/reference/ws_mad_logratio.md):
  median-centred log-ratio with optional per-sample MAD scaling; not
  Aitchison-valid but a robust descriptive transform for downstream
  learners that do not require zero-sum compositional coordinates.

- [`ws_dominance_score`](https://kstawiski.github.io/OmicSelector/reference/ws_dominance_score.md):
  panel-internal scalar QC dominance score reporting the fraction of
  total panel signal carried by the top-K features; flags pre-analytical
  failures (haemolysis, incomplete extraction) when calibrated against a
  Tier R reference distribution.

Helper utilities are provided alongside:
[`ws_default_sbp`](https://kstawiski.github.io/OmicSelector/reference/ws_default_sbp.md)
returns the curated v1 sequential-binary-partition tree for circulating
miRNA panels, and
[`ws_default_pivot_pool`](https://kstawiski.github.io/OmicSelector/reference/ws_default_pivot_pool.md)
returns the curated v1 ALR pivot pool (six low-variance miRNAs;
miR-92a-3p excluded because recent serum/plasma literature treats it as
erythrocyte-derived).

## Details

Background: circulating miRNA biomarkers were proposed as minimally
invasive cancer markers in 2008 (Mitchell et al., PNAS). In the
intervening years no clinically actionable cancer-screening test has
emerged from this source despite thousands of "candidate panel" papers.
The reasons are predominantly technical (batch effects, hemolysis
contamination, lack of standards, reference-cohort dependency,
cross-cohort leakage, absence of provenance/claim-gating) rather than
biological. These five within-sample compositional methods target two of
those failure modes directly: reference-cohort dependency (each method
scores a single sample using only its own panel values) and
contamination by dominant erythrocyte-derived species (miR-451a,
miR-486-5p, miR-16, miR-144-3p, miR-92a-3p - all explicitly handled by
either exclusion from the rCLR centering subset or by having dedicated
balances in the SBP tree).

Module A is one of six modules in the OmicSelector v2.4 toolkit. The
other modules - Module B (hemolysis correction), Module C (frozen batch
correction), Module D (outlier detection and conformal claim-gating),
Module E (multi-cancer / domain generalisation), and Module F
(self-supervised pretraining and ablation) - are documented in Paper 3
Methods and the per-module help pages.

## References

Aitchison J. (1986) The Statistical Analysis of Compositional Data.
Chapman & Hall.

Egozcue J. J., Pawlowsky-Glahn V. (2005) Groups of parts and their
balances in compositional data analysis. *Mathematical Geology* 37(7):
795-828.

Vandeputte D., Kathagen G., D'hoe K., et al. (2017) Quantitative
microbiome profiling links gut community variation to microbial load.
*Nature* 551(7681): 507-511. (Origin of the "detected_only" rCLR
convention.)

Mitchell P. S., Parkin R. K., Kroh E. M., et al. (2008) Circulating
microRNAs as stable blood-based markers for cancer detection.
*Proceedings of the National Academy of Sciences USA* 105(30):
10513-10518.

Stawiski K. (in preparation) Provenance-aware within-sample scoring for
circulating-microRNA biomarkers across cancers and platforms (Paper 3 of
the OmicSelector programme; Nature Methods target).
