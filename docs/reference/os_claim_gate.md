# Claim Gate for Panel-Null Benchmarks

Applies a fail-closed panel-claim policy: if the disease contrast has
not passed the provenance/identifiability review, no method can receive
a positive primary label. If identifiability passed, p-values are Holm
adjusted and a method passes only when adjusted p-value, null95 margin,
and null QC all satisfy the prespecified thresholds.

## Usage

``` r
os_claim_gate(
  results,
  identifiability_status,
  required_margin = 0.02,
  alpha = 0.05
)
```

## Arguments

- results:

  Data frame with columns `method_id`, `p_value`, `observed_stat`,
  `null_p95`, and `null_qc_pass`.

- identifiability_status:

  One of `"PASS"`, `"FAIL"`, or `"INCONCLUSIVE"`.

- required_margin:

  Numeric margin required over null95.

- alpha:

  Family-wise alpha for Holm-adjusted p-values.

## Value

Data frame with adjusted p-values, margins, and labels.
