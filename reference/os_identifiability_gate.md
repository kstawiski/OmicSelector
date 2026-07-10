# Apply a Fail-Closed Identifiability Gate

Adjudicates whether a public-omics contrast is identifiable enough to
permit downstream positive disease-signal language. Any missing or
failing criterion yields `"FAIL-DEMOTED"`.

## Usage

``` r
os_identifiability_gate(
  provenance_floor = NULL,
  provenance_auc = NULL,
  provenance_auc_ci_hi = NULL,
  single_covariate_max_auc = NULL,
  sensitivity_auc = NULL,
  within_blocks = NULL,
  u1_u2_overlap_fraction = NULL,
  u1_only_cases = NULL,
  u1_only_controls = NULL,
  physical_qc_pass = FALSE,
  anchor_leakage_pass = FALSE,
  floor_auc_ceiling = 0.65,
  floor_ci_ceiling = 0.7,
  single_covariate_ceiling = 0.7,
  sensitivity_auc_ceiling = 0.7,
  max_overlap_fraction = 0.8,
  min_u1_only_cases = 30L,
  min_u1_only_controls = 30L
)
```

## Arguments

- provenance_floor:

  Optional `os_provenance_floor` object.

- provenance_auc:

  Provenance-only AUC.

- provenance_auc_ci_hi:

  Upper confidence bound for provenance-only AUC.

- single_covariate_max_auc:

  Maximum single-provenance-covariate AUC.

- sensitivity_auc:

  Sensitivity-model AUC, for example a random-forest floor.

- within_blocks:

  Data frame from `os_within_provenance_blocks`.

- u1_u2_overlap_fraction:

  Fractional overlap of disease biological IDs.

- u1_only_cases:

  Count of disease biological IDs unique to U1.

- u1_only_controls:

  Count of control biological IDs unique to U1.

- physical_qc_pass:

  Logical physical-provenance and QC comparability flag.

- anchor_leakage_pass:

  Logical external-anchor leakage flag.

- floor_auc_ceiling:

  Maximum allowed provenance-floor AUC.

- floor_ci_ceiling:

  Maximum allowed upper CI for provenance-floor AUC.

- single_covariate_ceiling:

  Maximum allowed single-covariate AUC.

- sensitivity_auc_ceiling:

  Maximum allowed sensitivity-model AUC.

- max_overlap_fraction:

  Maximum allowed U1/U2 disease-ID overlap.

- min_u1_only_cases:

  Minimum U1-only cases.

- min_u1_only_controls:

  Minimum U1-only controls.

## Value

An `os_identifiability_gate` object.
