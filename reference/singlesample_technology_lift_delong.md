# Paired DeLong SE for technology-aware transfer lift

Computes the paired DeLong AUC-difference SE for a technology-aware
transfer method against its paired rCLR baseline on the same held-out
samples. This uses \`orient = "auc"\` to match the technology-aware
transfer convention \`max(AUC, 1 - AUC)\`.

## Usage

``` r
singlesample_technology_lift_delong(
  y,
  method_scores,
  baseline_scores,
  fold = NULL
)
```

## Arguments

- y:

  Binary held-out outcome vector.

- method_scores:

  Numeric method scores on the held-out samples.

- baseline_scores:

  Numeric paired baseline scores on the same samples.

- fold:

  Optional fold identifier for mean-of-folds aggregation.

## Value

A list from \`singlesample_paired_auc_diff_se()\`.
