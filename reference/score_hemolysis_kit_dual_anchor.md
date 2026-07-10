# Score a panel against the dual hemolysis plus kit-stable denominator.

Score a panel against the dual hemolysis plus kit-stable denominator.

## Usage

``` r
score_hemolysis_kit_dual_anchor(
  expr,
  meta,
  panel,
  dual_fit = NULL,
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

- dual_fit:

  Optional fit from \`fit_hemolysis_kit_dual_anchor()\`.

- train_expr:

  Optional training matrix used to fit \`dual_fit\`.

- train_meta:

  Optional training metadata used with \`train_expr\`.

- kit_col:

  Optional kit column.

- min_eval_anchors:

  Minimum anchors required in evaluation data.

- pseudocount:

  Optional additive pseudocount before log transform.

- ...:

  Additional arguments passed to \`fit_hemolysis_kit_dual_anchor()\`.
