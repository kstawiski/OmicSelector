# BenchmarkService: Nested Cross-Validation with Zero Leakage

R6 class that enforces proper nested cross-validation for biomarker
discovery. Implements the outer loop (evaluation) and inner loop
(selection) pattern required for unbiased performance estimation.

## Details

The BenchmarkService guarantees scientific validity by: - Enforcing that
feature selection occurs in the inner loop only - Computing the Nogueira
Stability Index across outer folds - Tracking which features are
selected in each fold for consensus analysis - Preventing any access to
test data during training/selection

## Methods

### Public methods

- [`BenchmarkService$new()`](#method-BenchmarkService-new)

- [`BenchmarkService$add_learner()`](#method-BenchmarkService-add_learner)

- [`BenchmarkService$run()`](#method-BenchmarkService-run)

- [`BenchmarkService$get_results()`](#method-BenchmarkService-get_results)

- [`BenchmarkService$print()`](#method-BenchmarkService-print)

- [`BenchmarkService$clone()`](#method-BenchmarkService-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new BenchmarkService

#### Usage

    BenchmarkService$new(
      task,
      outer_folds = 5,
      inner_folds = 3,
      stratify = TRUE,
      groups = NULL,
      seed = NULL
    )

#### Arguments

- `task`:

  An mlr3 Task or OmicPipeline object

- `outer_folds`:

  Number of outer CV folds (evaluation)

- `inner_folds`:

  Number of inner CV folds (selection/tuning)

- `stratify`:

  Logical, whether to stratify by outcome

- `groups`:

  Optional column name for grouped CV (e.g., patient_id)

- `seed`:

  Random seed for reproducibility

#### Returns

A new BenchmarkService object

------------------------------------------------------------------------

### Method `add_learner()`

Add a learner to benchmark

#### Usage

    BenchmarkService$add_learner(learner, id = NULL)

#### Arguments

- `learner`:

  A Learner, GraphLearner, or AutoFSelector

- `id`:

  Optional identifier for the learner

#### Returns

Self (for chaining)

------------------------------------------------------------------------

### Method `run()`

Run the nested cross-validation benchmark

#### Usage

    BenchmarkService$run(measures = NULL, parallel = FALSE)

#### Arguments

- `measures`:

  List of performance measures (default: AUC, accuracy)

- `parallel`:

  Logical, whether to run in parallel

#### Returns

A NestedCVResult object

------------------------------------------------------------------------

### Method `get_results()`

Get the most recent results

#### Usage

    BenchmarkService$get_results()

#### Returns

The NestedCVResult object

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print method

#### Usage

    BenchmarkService$print()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    BenchmarkService$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create benchmark service
service <- BenchmarkService$new(
  task = my_task,
  outer_folds = 5,
  inner_folds = 3
)

# Add learners with embedded feature selection
service$add_learner(my_graph_learner)

# Run nested CV
result <- service$run()

# Get stability metrics
stability <- result$get_stability()
} # }
```
