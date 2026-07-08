# Hurdle Filter for Zero-Inflated Data

Two-part filter designed for zero-inflated count/expression data. Tests
both the difference in zero-frequency AND the difference in magnitude
among non-zero values.

## Details

The Hurdle filter combines two tests: 1. \*\*Binary component\*\*:
Chi-square test on zero vs non-zero contingency table 2. \*\*Continuous
component\*\*: Wilcoxon test on non-zero values

The score is computed using Fisher's method to combine p-values:
`chi2 = -2 * sum(log(p_values))`

This captures genes that are good biomarkers because: - They are
detectable in one class but not the other (zero-frequency difference) -
They are upregulated in one class among expressed samples (magnitude
difference)

## Dictionary

This Filter can be instantiated via the dictionary \`mlr_filters\` or
with the associated sugar function \`flt()\`: “\` flt("hurdle")
mlr_filters\$get("hurdle") “\`

## Super class

[`mlr3filters::Filter`](https://mlr3filters.mlr-org.com/reference/Filter.html)
-\> `FilterHurdle`

## Methods

### Public methods

- [`FilterHurdle$new()`](#method-FilterHurdle-new)

- [`FilterHurdle$clone()`](#method-FilterHurdle-clone)

Inherited methods

- [`mlr3filters::Filter$calculate()`](https://mlr3filters.mlr-org.com/reference/Filter.html#method-calculate)
- [`mlr3filters::Filter$format()`](https://mlr3filters.mlr-org.com/reference/Filter.html#method-format)
- [`mlr3filters::Filter$help()`](https://mlr3filters.mlr-org.com/reference/Filter.html#method-help)
- [`mlr3filters::Filter$print()`](https://mlr3filters.mlr-org.com/reference/Filter.html#method-print)

------------------------------------------------------------------------

### Method `new()`

Creates a new instance of this Filter.

#### Usage

    FilterHurdle$new()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    FilterHurdle$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
