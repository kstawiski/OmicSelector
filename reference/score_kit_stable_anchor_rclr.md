# Score a panel against train-only kit-stable anchors.

Score a panel against train-only kit-stable anchors.

## Usage

``` r
score_kit_stable_anchor_rclr(
  expr,
  meta,
  panel,
  anchor_fit = NULL,
  train_expr = NULL,
  train_meta = NULL,
  kit_col = NULL,
  min_eval_anchors = 3L,
  pseudocount = NULL,
  ...
)
```

## Arguments

- expr:

  Numeric matrix, samples x features.

- meta:

  Data frame aligned to \`expr\`.

- panel:

  Character vector of panel features.

- anchor_fit:

  Optional fit from \`fit_kit_stable_anchors()\`.

- train_expr:

  Optional training matrix used to fit anchors when \`anchor_fit\` is
  NULL.

- train_meta:

  Optional training metadata used with \`train_expr\`.

- kit_col:

  Optional kit column.

- min_eval_anchors:

  Minimum anchors required in evaluation data.

- pseudocount:

  Optional additive pseudocount before log transform.

- ...:

  Additional arguments passed to \`fit_kit_stable_anchors()\`.
