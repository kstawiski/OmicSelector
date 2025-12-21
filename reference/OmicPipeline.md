# OmicPipeline: Zero-Leakage Feature Selection Pipeline

R6 class that encapsulates the mlr3 pipeline for biomarker discovery.
Guarantees zero data leakage by enforcing all preprocessing, feature
selection, and model training within proper cross-validation folds.

## Details

OmicPipeline is the central class for OmicSelector 2.0. It replaces the
legacy script-based approach with a rigorous, composable, and
reproducible architecture.

Key features: - All preprocessing (imputation, scaling) occurs inside CV
folds - Feature selection is embedded in the inner loop of nested CV -
Oversampling (SMOTE/ROSE) is applied only to training data per fold -
Factory methods generate configured GraphLearners

## Methods

### Public methods

- [`OmicPipeline$new()`](#method-OmicPipeline-new)

- [`OmicPipeline$create_graph_learner()`](#method-OmicPipeline-create_graph_learner)

- [`OmicPipeline$create_auto_fselector()`](#method-OmicPipeline-create_auto_fselector)

- [`OmicPipeline$benchmark()`](#method-OmicPipeline-benchmark)

- [`OmicPipeline$get_task()`](#method-OmicPipeline-get_task)

- [`OmicPipeline$get_feature_names()`](#method-OmicPipeline-get_feature_names)

- [`OmicPipeline$is_multi_omics()`](#method-OmicPipeline-is_multi_omics)

- [`OmicPipeline$get_modality_info()`](#method-OmicPipeline-get_modality_info)

- [`OmicPipeline$get_modality_features()`](#method-OmicPipeline-get_modality_features)

- [`OmicPipeline$print()`](#method-OmicPipeline-print)

- [`OmicPipeline$clone()`](#method-OmicPipeline-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new OmicPipeline object

#### Usage

    OmicPipeline$new(
      data,
      target,
      positive = NULL,
      patient_id = NULL,
      batch = NULL,
      id = "omic_task"
    )

#### Arguments

- `data`:

  Either a data.frame or a named list of data.frames for multi-omics.
  For multi-omics, use named list: list(rna = rna_data, mirna =
  mirna_data). Features will be namespaced: rna::gene1,
  mirna::hsa-miR-21.

- `target`:

  Name of the target column

- `positive`:

  Positive class label (for binary classification)

- `patient_id`:

  Optional column name for patient grouping (prevents leakage)

- `batch`:

  Optional column name for batch information

- `id`:

  Unique identifier for this pipeline

#### Returns

A new OmicPipeline object

------------------------------------------------------------------------

### Method `create_graph_learner()`

Create a GraphLearner with proper leakage prevention

#### Usage

    OmicPipeline$create_graph_learner(
      filter = "anova",
      model = "ranger",
      n_features = 20,
      impute_method = "median",
      scale = TRUE,
      oversample = NULL,
      batch_correct = FALSE
    )

#### Arguments

- `filter`:

  Filter method name (e.g., "anova", "mrmr", "correlation")

- `model`:

  Model type (e.g., "ranger", "glmnet", "svm")

- `n_features`:

  Number of features to select (or proportion if \< 1)

- `impute_method`:

  Imputation method ("median", "mean", "sample")

- `scale`:

  Logical, whether to scale features

- `oversample`:

  Oversampling method (NULL, "smote", "rose")

- `batch_correct`:

  Logical or character. If TRUE, adds FrozenComBat batch correction
  using the batch column specified in pipeline creation. If a character
  string, uses that as the batch column name. Default: FALSE.

#### Returns

A mlr3 GraphLearner object

------------------------------------------------------------------------

### Method `create_auto_fselector()`

Create an AutoFSelector for inner-loop feature selection tuning

#### Usage

    OmicPipeline$create_auto_fselector(
      learner,
      filter_values = c(5, 10, 20, 50),
      inner_resampling = NULL,
      measure = NULL
    )

#### Arguments

- `learner`:

  A Learner or GraphLearner

- `filter_values`:

  Vector of n_features values to try

- `inner_resampling`:

  Inner resampling strategy

- `measure`:

  Performance measure

#### Returns

An AutoFSelector object

------------------------------------------------------------------------

### Method `benchmark()`

Run benchmark with proper nested cross-validation

#### Usage

    OmicPipeline$benchmark(learners, outer_folds = 5, stratify = TRUE, seed = NULL)

#### Arguments

- `learners`:

  List of learners to benchmark

- `outer_folds`:

  Number of outer CV folds

- `stratify`:

  Logical, whether to stratify by outcome

- `seed`:

  Random seed for reproducibility

#### Returns

A BenchmarkResult object

------------------------------------------------------------------------

### Method `get_task()`

Get the underlying mlr3 Task

#### Usage

    OmicPipeline$get_task()

#### Returns

The mlr3 Task object

------------------------------------------------------------------------

### Method `get_feature_names()`

Get feature names

#### Usage

    OmicPipeline$get_feature_names()

#### Returns

Character vector of feature names (namespaced for multi-omics)

------------------------------------------------------------------------

### Method `is_multi_omics()`

Check if this is a multi-omics pipeline

#### Usage

    OmicPipeline$is_multi_omics()

#### Returns

Logical

------------------------------------------------------------------------

### Method [`get_modality_info()`](https://biostat.umed.pl/OmicSelector/reference/get_modality_info.md)

Get modality information for multi-omics data

#### Usage

    OmicPipeline$get_modality_info()

#### Returns

A data.frame with modality details, or NULL for single-modality

------------------------------------------------------------------------

### Method `get_modality_features()`

Get features for a specific modality

#### Usage

    OmicPipeline$get_modality_features(modality)

#### Arguments

- `modality`:

  Name of the modality (e.g., "rna", "mirna")

#### Returns

Character vector of feature names for that modality

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print method

#### Usage

    OmicPipeline$print()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    OmicPipeline$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create pipeline from data
pipeline <- OmicPipeline$new(
  data = my_data,
  target = "outcome",
  positive = "Case"
)

# Create a graph learner with feature selection
learner <- pipeline$create_graph_learner(
  filter = "anova",
  model = "ranger",
  n_features = 20
)

# Run nested cross-validation
result <- pipeline$benchmark(learner, outer_folds = 5, inner_folds = 3)
} # }
```
