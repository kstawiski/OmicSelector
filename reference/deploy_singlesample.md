# Build a frozen single-sample deployment scorer

Fits one rostered single-sample method and returns a frozen deployment
object that can score one incoming specimen without a scored batch,
co-resident reference cohort, or batch-correction step. This is a
deployability API, not a benchmark-ranking claim: the default
`"ws-balance-ilr"` is a reasonable compositional default, not a
certified winner. On the 21-cohort benchmark no single-sample method
robustly cleared +0.05 over trimmed-rCLR, and cross-platform transfer
was null.

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
  verify = TRUE
)
```

## Arguments

- X_train:

  Numeric matrix or data frame of training specimens (samples x
  features).

- y_train:

  Numeric/integer 0/1 labels aligned to `X_train`.

- method:

  Method identifier from
  [`singlesample_method_roster`](https://kstawiski.github.io/OmicSelector/reference/singlesample_method_roster.md).

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

## Value

An S3 object of class `singlesample_deployable` containing the frozen
`model`, `method_id`, resolved `score_fn`, `n_train`, `fit_time`, and
roster `meta`. Its value is deployability: scoring an incoming specimen
from that specimen plus frozen fit-time parameters only.

## See also

[`score_specimen`](https://kstawiski.github.io/OmicSelector/reference/score_specimen.md),
[`is_singlesample_deployable`](https://kstawiski.github.io/OmicSelector/reference/is_singlesample_deployable.md),
[`singlesample_method_roster`](https://kstawiski.github.io/OmicSelector/reference/singlesample_method_roster.md)
