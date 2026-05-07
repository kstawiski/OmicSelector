# OmicSelector_propensity_score_matching

Lightweight propensity score matching (nearest neighbor) using base R.

## Usage

``` r
OmicSelector_propensity_score_matching(
  dataset,
  match_by = c("age_at_diagnosis", "gender.x"),
  class_col = "Class",
  positive = "Case",
  method = "nearest",
  distance = "logit",
  ratio = 1,
  replace = FALSE,
  caliper = NULL,
  na_action = c("omit", "median")
)
```

## Arguments

- dataset:

  Original dataset.

- match_by:

  Character vector of covariate names to match on.

- class_col:

  Outcome column name (default: "Class").

- positive:

  Label for the treated/positive class (default: "Case").

- method:

  Matching method (currently only "nearest").

- distance:

  Propensity model type (currently only "logit").

- ratio:

  Number of controls per treated (default: 1).

- replace:

  Logical; allow matching controls with replacement.

- caliper:

  Optional maximum absolute propensity score distance.

- na_action:

  How to handle missing values in match_by: "omit" or "median".

## Value

Matched dataset (treated + matched controls).
