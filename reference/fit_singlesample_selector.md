# Fit a leakage-resistant single-sample method selector

Fits an experimental meta-scorer that selects or combines registered
single-sample methods using grouped cross-fitted predictions from the
supplied training data only. The fitted object scores each new specimen
from frozen base models; it never uses a co-resident test batch or test
outcomes.

Three prespecified routes are available:

- \`best_auc\`:

  Select the method with the largest mean grouped inner-fold AUC.

- \`lower_auc_ci\`:

  Select the method with the largest lower 95 percent stratified
  group-bootstrap limit for pooled cross-fitted AUC.

- \`simplex_stack\`:

  Fit a non-negative logistic stack whose weights sum to one, using
  cross-fitted group-level scores and a fixed numerical stabilizer. No
  stack hyperparameter is tuned.

This API establishes a leakage-resistant fitting contract; it is not a
claim that selection or stacking improves discrimination. Such a claim
requires a separate outer validation.

## Usage

``` r
fit_singlesample_selector(
  X_train,
  y_train,
  methods = NULL,
  route = c("best_auc", "lower_auc_ci", "simplex_stack", "all"),
  group_id = NULL,
  meta_train = NULL,
  annotation = NULL,
  inner_folds = 5L,
  seed = 1L,
  bootstrap_reps = 1000L,
  verify_base = TRUE,
  stack_stabilizer = 1e-06
)

# S3 method for class 'singlesample_selector_set'
print(x, ...)

# S3 method for class 'singlesample_selector_ineligible'
print(x, ...)
```

## Arguments

- X_train:

  Numeric matrix or data frame (specimens x features).

- y_train:

  Binary training outcome.

- methods:

  Character method identifiers. \`NULL\` uses
  \[singlesample_selector_candidates()\].

- route:

  One of \`"best_auc"\`, \`"lower_auc_ci"\`, \`"simplex_stack"\`, or
  \`"all"\`. The \`"all"\` route reuses one cross-fitting pass and
  returns a \`singlesample_selector_set\` containing all three
  prespecified routes.

- group_id:

  Biological/provenance group per row. Rows in one group are never split
  across inner folds. Defaults to one group per row.

- meta_train:

  Optional training metadata, one row per specimen.

- annotation:

  Optional feature annotation passed to compatible base methods.

- inner_folds:

  Requested number of grouped inner folds. It is reduced to the smaller
  class-specific group count and must remain at least two.

- seed:

  Integer seed for fold construction and deterministic derived
  per-method seeds.

- bootstrap_reps:

  Number of stratified biological-group bootstrap draws for
  \`lower_auc_ci\`.

- verify_base:

  Whether \[deploy_singlesample()\] should run its row- equivariance
  gate for every base fit.

- stack_stabilizer:

  Fixed non-negative L2 numerical stabilizer for the simplex stack. It
  is not tuned from outcomes.

- x:

  A \`singlesample_selector_set\` object.

- ...:

  Unused.

## Value

A \`singlesample_selector\`, \`singlesample_selector_set\`, or a
\`singlesample_selector_ineligible\` object with an explicit audit when
the requested route cannot be fitted.

## See also

\[score_singlesample_selector()\], \[deploy_singlesample()\]
