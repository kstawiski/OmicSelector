# Paired DeLong Non-Inferiority Test

Tests whether a candidate model's AUC is non-inferior to a reference
AUC, using a paired DeLong test on the same predictions.

\*\*Common bugs this prevents\*\*:

1.  Using the wrong reference AUC (e.g., superseded Phase 2.1 instead of
    current Phase 2.1b)

2.  Using an unpaired CI check instead of a formal paired DeLong test

3.  Comparing AUCs from different sample denominators

## Usage

``` r
test_noninferiority(
  response,
  predictor_candidate,
  predictor_reference,
  delta = 0.03,
  alpha = 0.05
)
```

## Arguments

- response:

  Binary outcome vector (0/1), same for both models

- predictor_candidate:

  Candidate model predictions

- predictor_reference:

  Reference model predictions. Must be on the SAME samples as
  predictor_candidate.

- delta:

  Non-inferiority margin (positive value, e.g., 0.03). Candidate is
  non-inferior if its AUC is within delta of the reference.

- alpha:

  Significance level (default 0.05 for one-sided test)

## Value

A list with:

- non_inferior:

  Logical: TRUE if candidate is non-inferior

- auc_candidate:

  Candidate AUC

- auc_reference:

  Reference AUC

- delta_auc:

  Candidate - Reference AUC

- delta:

  Non-inferiority margin used

- se_diff:

  Standard error of AUC difference (paired DeLong)

- z_stat:

  Z-statistic for non-inferiority

- p_value:

  One-sided p-value

- ci_lower_diff:

  Lower bound of 95% CI for AUC difference

## Examples

``` r
if (FALSE) { # \dontrun{
result <- test_noninferiority(
  response = y_test,
  predictor_candidate = pred_reduced_panel,
  predictor_reference = pred_full_panel,
  delta = 0.03
)
cat("Non-inferior:", result$non_inferior, "\n")
} # }
```
