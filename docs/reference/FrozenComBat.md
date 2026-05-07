# FrozenComBat R6 Class

R6 class for frozen batch effect correction using ComBat algorithm.
Parameters are estimated on training data and can be applied to new data
without re-estimation, preventing data leakage.

## Methods

### Public methods

- [`FrozenComBat$new()`](#method-FrozenComBat-new)

- [`FrozenComBat$fit()`](#method-FrozenComBat-fit)

- [`FrozenComBat$transform()`](#method-FrozenComBat-transform)

- [`FrozenComBat$fit_transform()`](#method-FrozenComBat-fit_transform)

- [`FrozenComBat$is_fitted()`](#method-FrozenComBat-is_fitted)

- [`FrozenComBat$get_batch_levels()`](#method-FrozenComBat-get_batch_levels)

- [`FrozenComBat$print()`](#method-FrozenComBat-print)

- [`FrozenComBat$clone()`](#method-FrozenComBat-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new FrozenComBat object

#### Usage

    FrozenComBat$new(parametric = TRUE, mean_only = FALSE)

#### Arguments

- `parametric`:

  Logical, use parametric priors (default TRUE)

- `mean_only`:

  Logical, only adjust means, not variances (default FALSE)

#### Returns

A FrozenComBat object

------------------------------------------------------------------------

### Method `fit()`

Fit ComBat parameters on training data

#### Usage

    FrozenComBat$fit(data, batch, covariates = NULL)

#### Arguments

- `data`:

  Numeric matrix with features in columns, samples in rows

- `batch`:

  Factor or character vector indicating batch membership

- `covariates`:

  Optional data.frame of biological covariates to preserve

#### Returns

Self (invisibly), for method chaining

------------------------------------------------------------------------

### Method [`transform()`](https://rdrr.io/r/base/transform.html)

Apply frozen ComBat correction to data

#### Usage

    FrozenComBat$transform(data, batch)

#### Arguments

- `data`:

  Numeric matrix with features in columns, samples in rows

- `batch`:

  Factor or character vector indicating batch membership

#### Returns

Batch-corrected matrix with same dimensions as input

------------------------------------------------------------------------

### Method `fit_transform()`

Fit and transform in one step (convenience method)

#### Usage

    FrozenComBat$fit_transform(data, batch, covariates = NULL)

#### Arguments

- `data`:

  Numeric matrix

- `batch`:

  Batch vector

- `covariates`:

  Optional covariates

#### Returns

Corrected matrix

------------------------------------------------------------------------

### Method `is_fitted()`

Check if the object has been fitted

#### Usage

    FrozenComBat$is_fitted()

#### Returns

Logical

------------------------------------------------------------------------

### Method `get_batch_levels()`

Get the batch levels learned during fitting

#### Usage

    FrozenComBat$get_batch_levels()

#### Returns

Character vector of batch level names

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print method Estimate ComBat parameters from training data This
implements the core empirical Bayes estimation from ComBat Apply frozen
correction parameters to new data

#### Usage

    FrozenComBat$print()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    FrozenComBat$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
