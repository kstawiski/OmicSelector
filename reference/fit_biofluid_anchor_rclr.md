# Fit biofluid-stable anchor rCLR

Selects denominator anchors with the lowest ratio of between-biofluid
variance to total variance in the training pool.

## Usage

``` r
fit_biofluid_anchor_rclr(
  train_expr,
  train_meta,
  panel,
  biofluid_col = "biofluid",
  pseudocount = NULL,
  anchor_n = 12L,
  min_anchor_n = 3L,
  min_biofluid_classes = 2L,
  exclude_features = .bfa_default_blood_cell_mirnas()
)
```

## Arguments

- train_expr:

  Numeric matrix, samples x features.

- train_meta:

  Data frame aligned to \`train_expr\`.

- panel:

  Character vector of panel features.

- biofluid_col:

  Column in \`train_meta\` with canonical or raw biofluid.

- pseudocount:

  Optional additive pseudocount before log transform.

- anchor_n:

  Number of anchor features to select.

- min_anchor_n:

  Minimum anchors required for an evaluable fit.

- min_biofluid_classes:

  Minimum biofluid classes in which an anchor must be detected.

- exclude_features:

  Features excluded from rCLR denominator selection.
