# Assert that a scorer is row-equivariant (single-sample deployable)

The decisive inclusion property for the manuscript method bank: at
inference a specimen's score must be computable from that specimen alone
given a frozen fitted model, with no dependence on which other specimens
are scored alongside it. This harness implements the strengthened
single-sample-deployability unit test prespecified in SAP Amendment \#4
(section 2): for the same specimen, the batch score must equal its
standalone score across ALL of

- \(a\) the singleton \\\\x_i\\\\ (the decisive single-sample case);

- \(b\) differently-composed batches (varied size/composition) and an
  adversarial batch (the specimen among extreme, scaled neighbours);

- \(c\) row-order permutations and batches with duplicated rows present;

- \(d\) no model-state mutation: the model's bytes are identical before
  and after scoring (catches BatchNorm running-stat updates / online
  recalibration).

Permutation alone cannot detect mean/quantile-centering coupling (column
means are permutation-invariant), which is why the singleton, subset,
duplicate, and adversarial cases always run.

## Usage

``` r
singlesample_assert_row_equivariant(
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

Invisibly \`TRUE\` on success; otherwise an error naming the failing
case.

## See also

\[singlesample_is_row_equivariant()\]
