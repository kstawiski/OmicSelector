# Hybrid Sequential Feature Selection (HSFS)

Implements a multi-stage feature reduction pipeline to handle
high-dimensional omics data efficiently. The pipeline chains: 1.
Variance thresholding (remove near-zero variance) 2. Univariate
filtering (ANOVA F-test) 3. RF importance filtering (single-pass
approximation of RFE) 4. LASSO regularization for final selection

Note: The RF importance stage is a single-pass approximation using
Random Forest variable importance, not true iterative Recursive Feature
Elimination.

This hybrid approach outperforms single-method selection for
high-dimensional data.

## References

Thelagathoti et al. "Hybrid feature selection approaches for machine
learning"

## Public fields

- `stages`:

  List of selection stages with parameters

- `selected_features`:

  Features selected at each stage

- `verbose`:

  Whether to print progress messages

## Methods

### Public methods

- [`SequentialSelector$new()`](#method-SequentialSelector-new)

- [`SequentialSelector$create_graph()`](#method-SequentialSelector-create_graph)

- [`SequentialSelector$create_learner()`](#method-SequentialSelector-create_learner)

- [`SequentialSelector$print()`](#method-SequentialSelector-print)

- [`SequentialSelector$clone()`](#method-SequentialSelector-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new SequentialSelector

#### Usage

    SequentialSelector$new(
      variance_threshold = 0.01,
      univariate_n = 5000,
      univariate_method = "anova",
      rfe_n = 1000,
      lasso_n = NULL,
      verbose = TRUE
    )

#### Arguments

- `variance_threshold`:

  Minimum variance to keep a feature (default: 0.01)

- `univariate_n`:

  Number of features after univariate filtering (default: 5000)

- `univariate_method`:

  Univariate method: "anova", "kruskal", "correlation"

- `rfe_n`:

  Number of features after RFE (default: 1000)

- `lasso_n`:

  Target number of features after LASSO (default: NULL = auto)

- `verbose`:

  Print progress messages

#### Returns

A new SequentialSelector object

------------------------------------------------------------------------

### Method `create_graph()`

Create an mlr3pipelines Graph for the sequential selection

#### Usage

    SequentialSelector$create_graph(task = NULL)

#### Arguments

- `task`:

  An mlr3 task to determine feature types

#### Returns

A Graph object that can be used in a GraphLearner

------------------------------------------------------------------------

### Method `create_learner()`

Create a complete GraphLearner with HSFS and a final classifier

#### Usage

    SequentialSelector$create_learner(
      model = "ranger",
      impute_method = "median",
      scale = TRUE
    )

#### Arguments

- `model`:

  Final classifier: "ranger", "xgboost", "glmnet", etc.

- `impute_method`:

  How to handle missing values: "median", "mean"

- `scale`:

  Whether to scale features

#### Returns

A GraphLearner ready for training/benchmarking

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print summary of the selector configuration

#### Usage

    SequentialSelector$print()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SequentialSelector$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
