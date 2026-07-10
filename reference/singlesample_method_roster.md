# Amendment \#4 method roster (frozen expansion)

Returns the frozen 74-row method roster for the OmicSelector
single-sample method-bank expansion (Statistical Analysis Plan Amendment
\#4, extended by Amendment \#5 which adds the Student-t-copula method
\`lrt-tcopula\` and Amendment \#6 which adds the KDE naive-Bayes method
\`lrt-nbkde\`). Unlike \[singlesample_method_bank()\], which is the
export-check view of functions already callable in the package, the
roster enumerates the full pre-specified set of single-sample scoring
methods (implemented and not-yet-implemented) together with the routing
metadata the benchmark instrument consumes.

## Usage

``` r
singlesample_method_roster()
```

## Value

A data frame with one row per frozen-roster method.

## Details

The roster is shipped as package data
(\`inst/extdata/singlesample_method_manifest.csv\`) so the R package and
the analysis instrument read one authoritative artifact. Each row
carries:

- method_id:

  Stable identifier (matches the analysis manifest).

- family:

  Human-readable grouping letter (A-O, NC).

- estimand:

  \`within\` or \`transfer\`; routes the SE/FDR family.

- role:

  \`baseline\`, \`discriminator\`, \`discriminator-conditional\`, or
  \`negative_control\`.

- tier:

  Execution substrate: \`R1\` (pure-R) or \`R2\` (reticulate).

- dep_route:

  Dependency / implementation route.

- fit_fn, score_fn:

  Package fit / score function names.

- pkg_status:

  \`present\`, \`primitive\`, \`new\`, or \`new-secondary\`.

- row_source:

  Instrument assembly path: \`per_cohort_rows\`,
  \`stratum_pooled_upstream\`, or \`transfer_fold_rows\`.

- lopbo_mechanism:

  \`within_cohort\`, \`transfer_fold\`, \`deferred\`, \`m4_decision\`
  (placeholder resolved at benchmark time), or \`na\`.

plus derived \`single_sample_gate\`, \`is_discriminator\`, and
\`is_negative_control\` columns.

## See also

\[singlesample_method_bank()\], \[singlesample_score_call()\],
\[singlesample_assert_row_equivariant()\]
