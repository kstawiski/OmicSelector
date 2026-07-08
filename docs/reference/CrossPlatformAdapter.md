# CrossPlatformAdapter R6 Class

An adapter for learning a source-platform representation and applying it
to target-platform data before model transfer.

## Value

A \`CrossPlatformAdapter\` object.

## Methods

### Public methods

- [`CrossPlatformAdapter$new()`](#method-CrossPlatformAdapter-new)

- [`CrossPlatformAdapter$fit()`](#method-CrossPlatformAdapter-fit)

- [`CrossPlatformAdapter$adapt()`](#method-CrossPlatformAdapter-adapt)

- [`CrossPlatformAdapter$get_source_data()`](#method-CrossPlatformAdapter-get_source_data)

- [`CrossPlatformAdapter$is_fitted()`](#method-CrossPlatformAdapter-is_fitted)

- [`CrossPlatformAdapter$print()`](#method-CrossPlatformAdapter-print)

- [`CrossPlatformAdapter$clone()`](#method-CrossPlatformAdapter-clone)

------------------------------------------------------------------------

### Method `new()`

Create a cross-platform adapter with the selected normalization strategy
and backend options.

#### Usage

    CrossPlatformAdapter$new(
      strategy = c("rank", "quantile", "zscore", "reference", "combat_pooled"),
      reference_features = NULL,
      reference_summary = c("mean", "median"),
      source_platform = "source",
      target_platform = "target",
      combat_parametric = TRUE,
      combat_mean_only = FALSE
    )

#### Arguments

- `strategy`:

  Adaptation strategy. One of \`"rank"\`, \`"quantile"\`, \`"zscore"\`,
  \`"reference"\`, or \`"combat_pooled"\`.

- `reference_features`:

  Optional character vector of reference miRNAs used when \`strategy =
  "reference"\`.

- `reference_summary`:

  Summary statistic for reference normalization: \`"mean"\` or
  \`"median"\`.

- `source_platform`:

  Label used for source samples in pooled ComBat.

- `target_platform`:

  Label used for target samples in pooled ComBat.

- `combat_parametric`:

  Logical, use parametric priors in pooled ComBat.

- `combat_mean_only`:

  Logical, only correct means in pooled ComBat.

------------------------------------------------------------------------

### Method `fit()`

Learn a source-platform reference representation.

#### Usage

    CrossPlatformAdapter$fit(source_data, source_labels = NULL)

#### Arguments

- `source_data`:

  Numeric matrix with samples in rows and features in columns.

- `source_labels`:

  Optional labels for the source data. Required for \`strategy =
  "combat_pooled"\`.

#### Returns

The adapter itself, invisibly.

------------------------------------------------------------------------

### Method `adapt()`

Adapt target-platform data to the learned source representation.

#### Usage

    CrossPlatformAdapter$adapt(target_data, target_labels = NULL)

#### Arguments

- `target_data`:

  Numeric matrix with samples in rows and features in columns.

- `target_labels`:

  Optional target labels. Required when \`strategy = "combat_pooled"\`.

#### Returns

Adapted target matrix with samples in rows and features in columns.

------------------------------------------------------------------------

### Method `get_source_data()`

Return the stored source matrix.

#### Usage

    CrossPlatformAdapter$get_source_data(adapted = TRUE)

#### Arguments

- `adapted`:

  Logical, return the adapted source matrix if \`TRUE\`, otherwise
  return the raw source matrix.

#### Returns

A numeric matrix.

------------------------------------------------------------------------

### Method `is_fitted()`

Check whether the adapter has been fitted.

#### Usage

    CrossPlatformAdapter$is_fitted()

#### Returns

Logical scalar.

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print a compact adapter summary.

#### Usage

    CrossPlatformAdapter$print()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    CrossPlatformAdapter$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
set.seed(42)
src <- matrix(rnorm(100 * 5), 100, 5)
tgt <- matrix(rnorm(50 * 5, mean = 2), 50, 5)

adapter <- CrossPlatformAdapter$new(strategy = "rank")
adapter$fit(src)
adapted_tgt <- adapter$adapt(tgt)
head(adapted_tgt)
#>      feature_1 feature_2 feature_3 feature_4 feature_5
#> [1,]      1.00      0.00      0.25      0.75      0.50
#> [2,]      1.00      0.25      0.75      0.00      0.50
#> [3,]      0.75      0.50      1.00      0.00      0.25
#> [4,]      0.25      1.00      0.75      0.50      0.00
#> [5,]      0.25      0.75      0.50      1.00      0.00
#> [6,]      0.50      0.75      1.00      0.00      0.25
```
