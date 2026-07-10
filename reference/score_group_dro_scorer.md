# Score held-out samples with a frozen Group-DRO scorer.

Score held-out samples with a frozen Group-DRO scorer.

## Usage

``` r
score_group_dro_scorer(fit, expr_matrix, sample_meta = NULL, ...)
```

## Arguments

- fit:

  Output from \`fit_group_dro_scorer()\`.

- expr_matrix:

  Numeric held-out matrix.

- sample_meta:

  Optional held-out metadata; used only for audit attributes, never for
  fitting.

- ...:

  Reserved for scorer-interface compatibility.

## Value

Numeric probability vector, higher meaning more disease-like.
