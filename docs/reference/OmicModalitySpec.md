# OmicModalitySpec R6 Class

Specification for a single omics modality including data, preprocessing,
and learner. Encapsulates everything needed to train/predict on one data
type.

## Public fields

- `id`:

  Unique identifier for this modality (e.g., "mRNA", "miRNA")

- `x`:

  Data matrix (samples x features)

- `learner`:

  mlr3 Learner object

- `preproc`:

  Optional mlr3pipelines Graph or PipeOp for preprocessing

## Methods

### Public methods

- [`OmicModalitySpec$new()`](#method-OmicModalitySpec-new)

- [`OmicModalitySpec$get_sample_ids()`](#method-OmicModalitySpec-get_sample_ids)

- [`OmicModalitySpec$n_features()`](#method-OmicModalitySpec-n_features)

- [`OmicModalitySpec$make_graph_learner()`](#method-OmicModalitySpec-make_graph_learner)

- [`OmicModalitySpec$as_task()`](#method-OmicModalitySpec-as_task)

- [`OmicModalitySpec$print()`](#method-OmicModalitySpec-print)

- [`OmicModalitySpec$clone()`](#method-OmicModalitySpec-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new modality specification

#### Usage

    OmicModalitySpec$new(id, x, learner, preproc = NULL)

#### Arguments

- `id`:

  Character, unique modality identifier

- `x`:

  Matrix or data.frame with samples in rows, features in columns.
  Rownames should be sample IDs.

- `learner`:

  mlr3 Learner or character shortcut (e.g., "classif.ranger")

- `preproc`:

  Optional preprocessing pipeline (Graph or PipeOp)

#### Returns

An OmicModalitySpec object

------------------------------------------------------------------------

### Method `get_sample_ids()`

Get sample IDs for this modality

#### Usage

    OmicModalitySpec$get_sample_ids()

#### Returns

Character vector of sample IDs

------------------------------------------------------------------------

### Method `n_features()`

Get number of features

#### Usage

    OmicModalitySpec$n_features()

#### Returns

Integer

------------------------------------------------------------------------

### Method `make_graph_learner()`

Create a GraphLearner combining preprocessing and learner

#### Usage

    OmicModalitySpec$make_graph_learner()

#### Returns

GraphLearner or Learner

------------------------------------------------------------------------

### Method `as_task()`

Create an mlr3 Task from this modality's data

#### Usage

    OmicModalitySpec$as_task(
      y,
      target_name = "target",
      task_type = c("classif", "regr"),
      sample_ids = NULL
    )

#### Arguments

- `y`:

  Target vector (must align with rownames of x)

- `target_name`:

  Name for target column

- `task_type`:

  "classif" or "regr"

- `sample_ids`:

  Optional subset of sample IDs to use

#### Returns

TaskClassif or TaskRegr

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print method

#### Usage

    OmicModalitySpec$print()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    OmicModalitySpec$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
