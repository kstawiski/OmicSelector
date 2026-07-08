# Plain-Text Bias-Audit Report

Convenience wrapper that runs
[`os_bias_audit`](https://kstawiski.github.io/OmicSelector/reference/os_bias_audit.md)
on a model/data bundle and returns the formatted report as a single
character string. Useful when the report needs to be captured into a
log, supplementary table, or pasted into a manuscript figure legend
without relying on the console print method.

## Usage

``` r
os_bias_audit_report(
  X = NULL,
  y = NULL,
  cohort = NULL,
  predictions = NULL,
  covariates = NULL,
  specimen_id = NULL,
  sample_id = NULL,
  feature_names = NULL,
  seed = 42L,
  audit = NULL
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

- audit:

  An optional pre-computed `os_bias_audit` object. If supplied, the
  audit is not recomputed and the remaining arguments are ignored. This
  makes it cheap to produce the text report alongside the structured
  result when both are needed.

## Value

A length-one character vector containing the full text report. The same
content that `print.os_bias_audit` would write to the console is
captured, with an additional attribute `"audit"` that holds the
underlying `os_bias_audit` object for downstream use.

## Examples

``` r
if (FALSE) { # \dontrun{
txt <- os_bias_audit_report(
  X = X, y = y, cohort = cohort,
  predictions = preds,
  covariates = list(age = age_years)
)
cat(txt)
writeLines(txt, "bias_audit.txt")
attr(txt, "audit")$bias_floor$apparent_auc
} # }
```
