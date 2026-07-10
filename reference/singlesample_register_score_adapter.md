# Register a non-canonical score adapter for a roster method

New roster methods follow the canonical convention \`score\_\*(model, X,
meta = NULL)\` and need no adapter. Legacy ("present") methods whose
score function uses a different argument order register a shim mapping
them onto the canonical \`function(model, X, meta)\` signature.

## Usage

``` r
singlesample_register_score_adapter(method_id, fn)
```

## Arguments

- method_id:

  Roster method identifier.

- fn:

  A function \`function(model, X, meta)\` returning one numeric score
  per row of \`X\`.

## Value

Invisibly \`TRUE\`.

## See also

\[singlesample_score_call()\]
