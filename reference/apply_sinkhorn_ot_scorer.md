# Apply a frozen Sinkhorn OT scorer.

Apply a frozen Sinkhorn OT scorer.

## Usage

``` r
apply_sinkhorn_ot_scorer(X_test, fit, meta_test = NULL)
```

## Arguments

- X_test:

  Samples x features matrix in the fitted feature space.

- fit:

  Object from \`fit_sinkhorn_ot_scorer()\`.

- meta_test:

  Optional metadata retained for a common scorer signature.

## Value

Projected CLR matrix with samples in rows and fit features in columns.
