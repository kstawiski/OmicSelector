# Per-Feature Batch-vs-Case Signal Partitioning

For each feature, computes three quantities that jointly distinguish
biological case signal from dataset/batch provenance:

- `F_dataset_healthy`: one-way ANOVA of feature expression against
  dataset identity, **using only healthy samples** (so case biology is
  removed);

- `F_case`: one-way ANOVA against case/control status, using all
  samples;

- `dataset_case_corr`: Pearson correlation between the sample-wise
  dataset-mean signature and the sample-wise case-mean signature.

Features with high `F_dataset_healthy` AND high `dataset_case_corr` are
the most suspect: their case signal is colinear with dataset identity.

## Usage

``` r
os_per_feature_batch_signal(X, y, cohort, feature_names = NULL)
```

## Arguments

- X:

  Numeric matrix or data frame (n samples x p features).

- y:

  Binary case vector (0/1). Length n.

- cohort:

  Factor or vector of cohort/dataset identity. Length n.

- feature_names:

  Optional character vector of feature labels.

## Value

A data frame with one row per feature, including a `suspect_score` that
combines healthy-only dataset ANOVA magnitude with `dataset_case_corr`,
ranked from most to least provenance-like.
