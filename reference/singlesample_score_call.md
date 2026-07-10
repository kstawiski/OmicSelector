# Score a roster method with the canonical single-sample interface

Resolves a roster method to its score function and invokes it through
the canonical adapter \`function(model, X_rows, meta_rows) -\>
numeric\`. New methods follow \`score\_\*(model, X, meta = NULL)\`
directly; legacy ("present") methods MUST be dispatched through a shim
registered with \[singlesample_register_score_adapter()\] (their
argument order predates the canonical convention, so calling them
positionally would mis-bind). Methods whose score function is not yet
implemented raise an explicit not-yet-implemented error.

## Usage

``` r
singlesample_score_call(
  method_id,
  model,
  X,
  meta = NULL,
  roster = singlesample_method_roster()
)
```

## Arguments

- method_id:

  Roster method identifier.

- model:

  A fitted model object produced by the method's \`fit\_\*\` function.

- X:

  Numeric matrix (samples x features), or a data frame coercible to one.

- meta:

  Optional per-sample metadata (data frame with \`nrow(X)\` rows).

- roster:

  The roster to resolve against; defaults to
  \[singlesample_method_roster()\].

## Value

Numeric vector of one score per row of \`X\`.

## See also

\[singlesample_method_roster()\],
\[singlesample_assert_row_equivariant()\]
