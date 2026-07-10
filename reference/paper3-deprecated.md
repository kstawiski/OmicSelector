# Deprecated \`paper3\_\*\` aliases

Back-compatibility aliases for the renamed \`singlesample\_\*\`
functions. Each forwards verbatim to its \`singlesample\_\*\`
counterpart and is retained so that code written against the former
\`paper3\_\*\` names keeps working. Prefer the \`singlesample\_\*\`
names; the \`paper3\_\*\` aliases may be removed in a future major
release.

## Usage

``` r
paper3_assert_method_bank_exports(...)

paper3_assert_row_equivariant(...)

paper3_bh_fdr_correct_blocked(...)

paper3_bh_fdr_correct_matched_null(...)

paper3_hanley_mcneil_auc_ci(...)

paper3_holm_correct_familywise(...)

paper3_is_row_equivariant(...)

paper3_make_loco_splits(...)

paper3_make_locto_splits(...)

paper3_matched_null_benchmark(...)

paper3_matched_null_benchmark_cv(...)

paper3_method_bank(...)

paper3_method_roster(...)

paper3_paired_auc_diff_se(...)

paper3_register_score_adapter(...)

paper3_score_call(...)

paper3_technology_lift_delong(...)
```

## Arguments

- ...:

  Arguments passed on to the corresponding \`singlesample\_\*\`
  function.

## Value

The value of the corresponding \`singlesample\_\*\` function.
