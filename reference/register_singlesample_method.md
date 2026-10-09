# Register a custom single-specimen scorer for this R session

Adds a fit/score pair to \[singlesample_method_roster()\],
\[deploy_singlesample()\] and \[singlesample_score_call()\].
Registration is local to the loaded package in this R session. Built-in
routes cannot be replaced. Deployment objects retain their resolved
custom scorer and strict support flag after a route is replaced or
removed.

## Usage

``` r
register_singlesample_method(
  method_id,
  fit_fn,
  score_fn,
  family = "custom",
  role = "discriminator",
  estimand = "custom single-specimen score",
  strict_support = FALSE,
  audit = TRUE,
  probe = NULL,
  overwrite = FALSE
)

unregister_singlesample_method(method_id)

registered_singlesample_methods()
```

## Arguments

- method_id:

  One syntactically valid R identifier, unused by built-in routes.
  Existing user routes require \`overwrite = TRUE\`.

- fit_fn:

  Function declaring \`X_train\` and \`y_train\`. Additional arguments
  must have defaults, except optional deployment inputs \`meta_train\`,
  \`meta\` and \`annotation\`. \`...\` is allowed but not required.

- score_fn:

  Function declaring \`model\`, \`X\` and \`meta = NULL\`, returning one
  finite numeric score per specimen. Additional arguments must have
  defaults; \`...\` is allowed. All fitted parameters must be stored in
  \`model\`. Scoring must not change the model or depend on other
  specimens or mutable external state.

- family, role, estimand:

  Non-empty character scalars describing the route. Use \`estimand =
  "within"\` to include a route in within-cohort benchmark helpers, or
  \`"transfer"\` for transfer routing. The descriptive default is
  excluded from within-cohort comparisons. Selector eligibility also
  requires role \`"baseline"\` or \`"discriminator"\`.

- strict_support:

  Logical scalar. If \`TRUE\`, \`score_fn\` must accept \`strict\`
  explicitly or through \`...\`; scoring passes this flag in both modes.
  The callback defines its strict feature policy. Otherwise strict
  requests fail with \`STRICT_UNSUPPORTED\`.

- audit:

  Logical scalar. If \`TRUE\`, fit the probe and run
  \[singlesample_assert_row_equivariant()\], including its serialized
  fitted- state immutability check. Both scoring modes are audited when
  strict support is declared. \`FALSE\` skips only registration
  auditing; deployment verification still defaults to \`TRUE\`.

- probe:

  Optional list with elements \`X\` (numeric matrix with unique feature
  names and at least four rows) and \`y\` (aligned binary 0/1 labels,
  both classes present), named or supplied in that order. The default
  has 12 positive-valued specimens and six features named \`feature1\`
  through \`feature6\`. Supply a probe for callbacks requiring a
  specific panel. Extra fit settings can be captured in a function
  wrapper with fixed configuration.

- overwrite:

  Logical scalar. Replace an existing user route only after all
  requested checks pass.

## Value

\`register_singlesample_method()\` returns invisibly \`TRUE\`.
\`unregister_singlesample_method()\` returns invisibly whether a user
route was removed. \`registered_singlesample_methods()\` returns a data
frame in identifier order with roster columns plus \`strict_support\`
and \`audited\`.

## Details

Registration failures are \`singlesample_input_error\` conditions with
codes \`REGISTRATION_INVALID_ID\`, \`REGISTRATION_ID_EXISTS\`,
\`REGISTRATION_INVALID_FUNCTION\`, \`REGISTRATION_INVALID_ARGUMENT\`, or
\`REGISTRATION_AUDIT_FAILED\`. Audit errors retain the original message.
Custom scoring uses the deployment input guards. Callback errors use
\`CUSTOM_SCORER_FAILED\`, unless already coded; invalid or non-finite
score output uses \`CUSTOM_INVALID_SCORES\`.

Auditing tests the supplied probe with absolute tolerance 1e-8. It is
not a proof for every possible input or evidence of clinical
performance. Models with external pointers cannot pass the default
serialization audit. Audit fitting preserves the caller's random seed,
but callbacks must avoid other external side effects. Re-register routes
in each new R session to fit new deployments; saved deployments contain
their scorer but still require its dependencies and a self-contained
callback environment.

## See also

\[deploy_singlesample()\], \[score_specimen()\]
