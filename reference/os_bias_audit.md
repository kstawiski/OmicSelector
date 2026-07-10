# One-Shot Bias Audit Report

Orchestrator that runs all bias-audit diagnostics on a model + data
bundle and returns a structured report. Intended to be emitted as a
first-class output of every OmicSelector biomarker analysis.

## Usage

``` r
os_bias_audit(
  X,
  y,
  cohort,
  predictions = NULL,
  covariates = NULL,
  specimen_id = NULL,
  sample_id = NULL,
  feature_names = NULL,
  seed = 42L
)
```

## Arguments

- X:

  Feature matrix (n x p).

- y:

  Binary outcome.

- cohort:

  Cohort/dataset identity per sample.

- predictions:

  Optional numeric vector of model predictions (n).

- covariates:

  Optional named list of numeric or factor covariates (e.g., list(age =
  age_years, sex = sex_factor)).

- specimen_id:

  Optional specimen-level IDs for duplicate detection.

- sample_id:

  Optional sample-level IDs used when reporting duplicate pairs. If
  omitted, row names of `X` are used when available.

- feature_names:

  Optional feature labels.

- seed:

  Random seed.

## Value

An `os_bias_audit` object (list) with fields: `bias_floor`,
`covariate_aucs`, `feature_signal`, `duplicates`, `provenance_summary`,
`classifier_vs_floor` (if predictions provided).
