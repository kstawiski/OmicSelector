# Test whether a scorer is row-equivariant (non-throwing)

Convenience wrapper around \[singlesample_assert_row_equivariant()\]
returning a logical instead of erroring. Used both for affirmative
checks and for the negative test that a deliberately batch-coupled
scorer (e.g. the batchwise Sinkhorn plan) is \*not\* row-equivariant.

## Usage

``` r
singlesample_is_row_equivariant(
  score_fun,
  model,
  X,
  meta = NULL,
  n_perm = 4L,
  tol = 1e-08,
  n_singletons = 4L,
  model_digest = NULL,
  seed = 1L
)
```

## Arguments

- score_fun:

  A function \`function(model, X, meta)\` returning one numeric score
  per row of \`X\`.

- model:

  A fitted model object.

- X:

  Numeric matrix (samples x features), at least 4 rows.

- meta:

  Optional per-sample metadata (data frame with \`nrow(X)\` rows).

- n_perm:

  Number of random row permutations to test.

- tol:

  Absolute tolerance for score agreement (default 1e-8).

- n_singletons:

  Number of distinct rows to test as singletons (always includes the
  first and last row). Pass \`nrow(X)\` for an exhaustive sweep.

- model_digest:

  Optional \`function(model)\` returning a comparable snapshot of the
  model's mutable state. Defaults to a serialization snapshot; supply a
  custom closure for torch / external-pointer models that do not
  serialize deterministically.

- seed:

  Random seed for reproducible permutations/subsets.

## Value

\`TRUE\` if row-equivariant, \`FALSE\` otherwise.

## See also

\[singlesample_assert_row_equivariant()\]
