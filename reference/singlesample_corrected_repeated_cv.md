# Corrected repeated-CV inference for a paired method effect

Summarises paired effects observed in the same repeated cross-validation
strata using the Nadeau–Bengio correction. The standard error is
\`sqrt((1 / m + mean(n_test / n_train)) \* var(effect))\`, with \`m -
1\` degrees of freedom. The function fails closed unless the observed
stratum identifiers match the complete prespecified set exactly.

## Usage

``` r
singlesample_corrected_repeated_cv(
  effect,
  n_test,
  n_train,
  stratum_id,
  expected_m,
  expected_strata,
  margin = 0.05,
  conf_level = 0.95
)
```

## Arguments

- effect:

  Numeric paired effect for every shared CV stratum, ordered as method A
  minus method B.

- n_test, n_train:

  Positive test and training counts for each stratum. Use
  biological-group counts for the group-collapsed primary estimand;
  profile counts are appropriate only for an explicit profile-weighted
  sensitivity analysis.

- stratum_id:

  Unique observed stratum identifiers.

- expected_m:

  Prespecified number of shared strata.

- expected_strata:

  Complete prespecified stratum identifiers. Their order may differ from
  \`stratum_id\`, but membership must match exactly.

- margin:

  Positive relevance margin. The default is \`0.05\` AUC units.

- conf_level:

  Confidence level for the nominal t interval.

## Value

A named list containing the estimate, corrected uncertainty, confidence
interval, zero and shifted-null p-values, correction factor, and exact
stratum accounting.

## Details

The returned directional tests are fixed before outcomes are inspected:
\`p_above_margin\` tests \`H0: effect \<= margin\`, while
\`p_below_minus_margin\` tests \`H0: effect \>= -margin\`.
\`p_zero_two_sided\` tests a zero mean effect. When the corrected
standard error is zero, the confidence interval remains the point
estimate but inferential p-values are \`NA\`, rather than treating a
degenerate repeated-CV sample as certain.

## References

Nadeau C, Bengio Y. Inference for the generalization error. Machine
Learning. 2003;52:239–281.
