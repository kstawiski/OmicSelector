# Fit the structurally-informed ILR single-sample scorer (ss-struct-ilr)

Builds a genomic-cluster / GC-tertile Sequential Binary Partition from
the training feature set + a miRNA `annotation` table, computes ILR
balances per specimen via
[`ws_balance_ilr`](https://kstawiski.github.io/OmicSelector/reference/ws_balance_ilr.md),
standardises them with train-only mean/sd, and learns a frozen
diagonal-LDA linear rule (per-balance standardised case-minus-control
mean difference). The returned object serialises every artefact needed
to score one new specimen with no test-batch information.

## Usage

``` r
fit_ss_struct_ilr(X_train, y_train, meta_train = NULL, annotation, hp = list())
```

## Arguments

- X_train:

  Numeric matrix (samples x features) of non-negative abundances;
  columns are uniquely named miRNA ids.

- y_train:

  0/1 labels aligned to `X_train` (1 = case).

- meta_train:

  Optional per-sample metadata; accepted for interface uniformity and
  ignored.

- annotation:

  data.frame with columns `feature` (miRNA id), `cluster`
  (genomic-cluster id), `gc` (GC fraction).

- hp:

  List of frozen hyperparameters: `pseudocount` (default `NULL` = the
  scale-invariant per-sample default of
  [`ws_balance_ilr`](https://kstawiski.github.io/OmicSelector/reference/ws_balance_ilr.md)),
  `aggregate` ("gmean"/"trimmed_gmean", default "gmean"),
  `min_balance_coverage` (default 0.5), `min_part` (default 1L).

## Value

A `ss_struct_ilr_model` list: frozen `sbp`, `balance_names`, `center`,
`scale`, `weights`, `feature_universe`, the resolved `hp`, and an
`annotation_features` record.

## See also

[`score_ss_struct_ilr`](https://kstawiski.github.io/OmicSelector/reference/score_ss_struct_ilr.md),
[`ws_balance_ilr`](https://kstawiski.github.io/OmicSelector/reference/ws_balance_ilr.md)
