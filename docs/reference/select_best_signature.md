# Select Best Biomarker Signature from Nested CV Results

Multi-objective signature selection with three modes: constrained 1SE
rule, weighted scoring, or Pareto frontier. Designed for clinical
biomarker panels where parsimony and stability are as important as raw
performance.

## Usage

``` r
select_best_signature(
  nested_result,
  mode = c("constrained_1se", "weighted", "pareto"),
  metric = "classif.auc",
  metric_higher_better = TRUE,
  auc_min = NULL,
  stability_min = NULL,
  k_max = NULL,
  weights = c(performance = 0.5, stability = 0.3, parsimony = 0.2),
  na_stability = c("exclude_if_constrained", "worst", "neutral")
)
```

## Arguments

- nested_result:

  A NestedCVResult object from BenchmarkService\$run()

- mode:

  Selection mode: "constrained_1se" (default), "weighted", or "pareto"

- metric:

  Character, performance metric column name (default: "classif.auc")

- metric_higher_better:

  Logical, TRUE if higher metric is better (default: TRUE)

- auc_min:

  Numeric, minimum acceptable AUC (hard constraint)

- stability_min:

  Numeric, minimum acceptable stability index (hard constraint)

- k_max:

  Integer, maximum acceptable feature count (hard constraint for
  clinical panels)

- weights:

  Named numeric vector for weighted mode: c(performance=, stability=,
  parsimony=)

- na_stability:

  How to handle NA stability: "exclude_if_constrained" (default),
  "worst" (treat as worst), or "neutral" (treat as middle)

## Value

A data.table with candidate signatures ranked by the selection criteria.
Contains columns: learner_id, mean_metric, se_metric, stability, mean_k,
score (for weighted), pareto (for pareto mode), selected (best candidate
flag)

## Examples

``` r
if (FALSE) { # \dontrun{
# Run nested CV
result <- benchmark_service$run()

# Select best with 1SE rule (conservative)
best <- select_best_signature(result, mode = "constrained_1se")

# Select with explicit weights
best <- select_best_signature(result, mode = "weighted",
  weights = c(performance = 0.5, stability = 0.3, parsimony = 0.2))

# Get Pareto frontier for manual inspection
pareto_set <- select_best_signature(result, mode = "pareto")
} # }
```
