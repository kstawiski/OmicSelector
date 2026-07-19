# Matched AUC difference across repeated CV strata

Computes fixed-direction AUC for two methods on exactly aligned held-out
rows and applies
[`singlesample_corrected_repeated_cv()`](https://kstawiski.github.io/OmicSelector/reference/singlesample_corrected_repeated_cv.md)
to the stratum-level paired differences. The default preserves the
historical group-collapsed estimand; `analysis_level = "profile"` gives
every held-out profile equal weight while retaining
provenance-group-safe folds. Larger scores must already indicate the
positive class according to a training-frozen orientation.

## Usage

``` r
singlesample_matched_pair_auc(
  y,
  score_a,
  score_b,
  stratum_id,
  sample_id,
  group_id,
  expected_strata,
  analysis_level = c("group", "profile"),
  margin = 0.05,
  conf_level = 0.95
)
```

## Arguments

- y:

  Binary held-out outcome.

- score_a, score_b:

  Finite held-out scores for methods A and B on the same rows.

- stratum_id:

  Repeated-CV stratum per row, such as \`seed::fold\`.

- sample_id:

  Held-out profile identifier. Each profile may occur only once in a
  stratum.

- group_id:

  Biological/provenance group identifier. Group labels must be
  homogeneous.

- expected_strata:

  Complete prespecified repeated-CV strata.

- analysis_level:

  Either `"group"` (collapse repeated profiles within
  biological/provenance group) or `"profile"` (equal profile weight).

- margin:

  Positive relevance margin in AUC units.

- conf_level:

  Confidence level for the nominal corrected t interval.

## Value

A list with one \`strata\` data frame and one corrected \`summary\`
list.
