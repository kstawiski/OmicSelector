# Export Model as Vetiver for Deployment

Wraps an mlr3 Learner or GraphLearner in a vetiver model object for
versioning and deployment. Vetiver provides model versioning, API
generation, and monitoring capabilities.

## Usage

``` r
export_vetiver(
  learner,
  model_name,
  task,
  description = NULL,
  metadata = NULL,
  board = NULL
)
```

## Arguments

- learner:

  A trained mlr3 Learner or GraphLearner

- model_name:

  Name for the vetiver model (required)

- task:

  The mlr3 Task used for training (for metadata extraction)

- description:

  Optional description string

- metadata:

  Optional list of additional metadata

- board:

  Optional pins board for immediate deployment. If NULL, returns the
  vetiver model without pinning.

## Value

A vetiver_model object (or invisibly if board is provided)

## Examples

``` r
if (FALSE) { # \dontrun{
library(vetiver)
library(pins)

# Train model
task <- tsk("iris")
learner <- lrn("classif.ranger", predict_type = "prob")
learner$train(task)

# Export as vetiver
v <- export_vetiver(
  learner = learner,
  model_name = "iris_classifier",
  task = task,
  description = "Random Forest classifier for iris dataset"
)

# Pin to a board
board <- board_folder("models", versioned = TRUE)
vetiver_pin_write(board, v)

# Or export and pin in one step
export_vetiver(learner, "iris_classifier", task, board = board)
} # }
```
