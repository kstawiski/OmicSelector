# Score a panel against a RIN-weighted frozen training reference profile.

Score a panel against a RIN-weighted frozen training reference profile.

## Usage

``` r
score_rin_weighted_score(
  expr,
  meta,
  panel,
  rin_fit = NULL,
  train_expr = NULL,
  train_meta = NULL,
  rin_col = NULL,
  pseudocount = NULL
)
```

## Arguments

- expr:

  Numeric matrix, samples x features.

- meta:

  Data frame aligned to \`expr\`.

- panel:

  Character vector of panel features.

- rin_fit:

  Optional fit from \`fit_rin_weighted_reference()\`.

- train_expr:

  Optional training matrix used to fit \`rin_fit\`.

- train_meta:

  Optional training metadata used with \`train_expr\`.

- rin_col:

  Optional RNA-integrity-number column.

- pseudocount:

  Optional additive pseudocount before log transform.
