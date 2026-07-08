# Frozen Plate-Median Correction

R6 class for fold-frozen per-plate (or per-batch) location-shift
correction. Subtracts the per-plate median per feature, then adds back
the grand median of training-plate medians, so each plate is centered
without changing the overall mean. Optionally divides by a per-plate
scale (IQR or MAD) before re-multiplying by the grand scale.

## Details

This is intentionally simpler than
[`FrozenComBat`](https://kstawiski.github.io/OmicSelector/reference/FrozenComBat.md):
it assumes a pure location-shift (or location-and-scale) batch model and
is single-plate deployable. For ddPCR-style data where plate effects are
largely additive in log-space this is often the right tool.

Use
[`FrozenComBat`](https://kstawiski.github.io/OmicSelector/reference/FrozenComBat.md)
when batch effects are believed to interact with biological covariates
and an empirical-Bayes shrinkage is desirable.

## Public methods

- `$new(scale)`:

  Constructor; `scale` ∈ {"none", "iqr", "mad"}.

- `$fit(data, plate)`:

  Learn per-plate medians (and optional scales) from training data only.

- `$transform(data, plate)`:

  Apply frozen parameters to new data. Plates not seen during fit fall
  back to grand-median centering and emit a warning.

- `$fit_transform(data, plate)`:

  Convenience for fit + transform.

- `$is_fitted()`:

  Logical.

- `$get_plate_levels()`:

  Character vector of training plate labels.

## Methods

### Public methods

- [`os_plate_median_frozen$new()`](#method-os_plate_median_frozen-new)

- [`os_plate_median_frozen$fit()`](#method-os_plate_median_frozen-fit)

- [`os_plate_median_frozen$transform()`](#method-os_plate_median_frozen-transform)

- [`os_plate_median_frozen$fit_transform()`](#method-os_plate_median_frozen-fit_transform)

- [`os_plate_median_frozen$is_fitted()`](#method-os_plate_median_frozen-is_fitted)

- [`os_plate_median_frozen$get_plate_levels()`](#method-os_plate_median_frozen-get_plate_levels)

- [`os_plate_median_frozen$clone()`](#method-os_plate_median_frozen-clone)

------------------------------------------------------------------------

### Method `new()`

Construct a new os_plate_median_frozen object.

#### Usage

    os_plate_median_frozen$new(scale = c("none", "iqr", "mad"))

#### Arguments

- `scale`:

  One of `"none"` (default), `"iqr"`, or `"mad"`.

------------------------------------------------------------------------

### Method `fit()`

Fit per-plate centering parameters on training data.

#### Usage

    os_plate_median_frozen$fit(data, plate)

#### Arguments

- `data`:

  Numeric matrix; samples in rows, features in columns.

- `plate`:

  Character or factor vector; length `nrow(data)`.

------------------------------------------------------------------------

### Method [`transform()`](https://rdrr.io/r/base/transform.html)

Apply frozen correction to new data.

#### Usage

    os_plate_median_frozen$transform(data, plate)

#### Arguments

- `data`:

  Numeric matrix.

- `plate`:

  Character or factor vector.

------------------------------------------------------------------------

### Method `fit_transform()`

Fit and transform in one call.

#### Usage

    os_plate_median_frozen$fit_transform(data, plate)

#### Arguments

- `data`:

  Numeric matrix.

- `plate`:

  Plate vector.

------------------------------------------------------------------------

### Method `is_fitted()`

Logical: has the object been fitted?

#### Usage

    os_plate_median_frozen$is_fitted()

------------------------------------------------------------------------

### Method `get_plate_levels()`

Character vector of plate labels seen during fit.

#### Usage

    os_plate_median_frozen$get_plate_levels()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    os_plate_median_frozen$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
set.seed(7)
plate <- rep(c("P1", "P2", "P3"), each = 10)
x <- matrix(rnorm(60, mean = 0, sd = 0.5), nrow = 30)
x[plate == "P2", ] <- x[plate == "P2", ] + 3
x[plate == "P3", ] <- x[plate == "P3", ] - 2
fc <- os_plate_median_frozen$new()
fc$fit(x[1:24, ], plate[1:24])
corrected <- fc$transform(x[25:30, ], plate[25:30])
```
