# Build a frozen single-sample deployment scorer

Fits one rostered single-sample method and returns a frozen deployment
object that can score one incoming specimen without a scored batch,
co-resident reference cohort, or batch-correction step. This is a
deployability API, not a benchmark-ranking claim: the default
`"ws-balance-ilr"` is a reasonable compositional default, not a
certified winner. It does not assert that any single-sample method
out-discriminates a batch-corrected pipeline; validate the method choice
on your own data.

Reticulate methods (roster tier R2; require a Python backend / venv) are
intentionally not fitted by this base-R deployment wrapper.

## Usage

``` r
deploy_singlesample(
  X_train,
  y_train,
  method = "ws-balance-ilr",
  meta_train = NULL,
  annotation = NULL,
  verify = TRUE,
  nondetect_limit = 0,
  preprocessor = NULL,
  group_train = NULL
)
```

## Arguments

- X_train:

  Numeric matrix or data frame of training specimens (samples x
  features). With a preprocessor, supply raw measurements.

- y_train:

  Numeric/integer 0/1 labels aligned to `X_train`.

- method:

  Method identifier from
  [`singlesample_method_roster`](https://kstawiski.github.io/OmicSelector/reference/singlesample_method_roster.md).
  Custom routes can be added with \[register_singlesample_method()\].

- meta_train:

  Optional training metadata, one row per training specimen.

- annotation:

  Optional feature annotation passed only to methods whose fit signature
  includes an `annotation` argument.

- verify:

  Logical; if `TRUE`, run
  [`singlesample_assert_row_equivariant`](https://kstawiski.github.io/OmicSelector/reference/singlesample_assert_row_equivariant.md)
  on the training rows (or a deterministic subset for large inputs) so
  construction fails loudly when the fitted scorer is not single-sample
  deployable.

- nondetect_limit:

  Finite non-negative abundance threshold, frozen at deployment
  construction (default 0). For an assay mapping undetected Ct 40 to
  abundance 1, use 1. This threshold is on the abundance scale after any
  preprocessing. An attached preprocessor also supplies its detection
  mask.

- preprocessor:

  Optional frozen `singlesample_preprocessor` from
  [`fit_singlesample_preprocessor`](https://kstawiski.github.io/OmicSelector/reference/fit_singlesample_preprocessor.md).
  The default `NULL` uses the existing abundance-input path. When
  supplied, training, scoring and behavioural verification all accept
  raw measurements.

- group_train:

  Optional biological-group identifiers, one per training row, without
  missing values. Whole groups stay together in inner resampling.
  \`NULL\` preserves the original row-based fit. Invalid or infeasible
  groups raise \`singlesample_input_error\` with code
  \`INVALID_GROUPS\`.

## Value

An S3 object of class `singlesample_deployable` containing the frozen
`model`, `method_id`, resolved `score_fn`, `n_train`, `fit_time`, roster
`meta`, and the `preprocessor` when supplied. Its value is
deployability: scoring an incoming specimen from that specimen plus
frozen fit-time parameters only.

## See also

[`score_specimen`](https://kstawiski.github.io/OmicSelector/reference/score_specimen.md),
[`is_singlesample_deployable`](https://kstawiski.github.io/OmicSelector/reference/is_singlesample_deployable.md),
[`singlesample_method_roster`](https://kstawiski.github.io/OmicSelector/reference/singlesample_method_roster.md)
