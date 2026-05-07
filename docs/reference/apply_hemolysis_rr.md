# Apply frozen hemolysis-RR correction to new samples

Subtracts the predicted hemolysis (and optional platelet) nuisance
contribution from one or more test samples using the frozen per-feature
regression coefficients from
[`fit_hemolysis_rr`](https://kstawiski.github.io/OmicSelector/reference/fit_hemolysis_rr.md).

Corrections are clipped to \\\pm\\`clip_max` log-units to prevent
over-correction. If `report_oob = TRUE` (default), an `oob` attribute
flags samples whose hemolysis / platelet scores fall outside the
training range.

## Usage

``` r
apply_hemolysis_rr(
  x,
  fit,
  hemo_score,
  platelet_score = NULL,
  report_oob = TRUE
)
```

## Arguments

- x:

  Numeric vector (single sample) or matrix (samples \\\times\\
  features).

- fit:

  Object of class `hemolysis_rr_fit` from
  [`fit_hemolysis_rr`](https://kstawiski.github.io/OmicSelector/reference/fit_hemolysis_rr.md).

- hemo_score:

  Numeric vector (length = number of samples).

- platelet_score:

  Optional numeric vector (same length as `hemo_score`).

- report_oob:

  Logical. If `TRUE`, attaches an `oob` attribute flagging
  out-of-training-range samples. Default `TRUE`.

## Value

Same shape as `x`. Hemolysis-corrected log-expression values, with
attributes `oob` (if `report_oob`) and `fit_summary`.
