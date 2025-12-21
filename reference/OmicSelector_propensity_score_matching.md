# OmicSelector_propensity_score_matching

Propensity score matching.

## Usage

``` r
OmicSelector_propensity_score_matching(
  dataset,
  match_by = c("age_at_diagnosis", "gender.x"),
  method = "nearest",
  distance = "logit"
)
```

## Arguments

- dataset:

  Original dataset.

- match_by:

  Vector describing by which variables should the dataset be matched.

- method:

  Passed to \`matchit()\`.

- distance:

  Passed to \`matchit()\`.
