# Zero-Proportion Filter

Simple filter based on the difference in zero-proportion between
classes. Useful for quickly identifying on/off switches in gene
expression.

## Details

Computes: \|P(X=0\|Class1) - P(X=0\|Class2)\| Higher values indicate
features that are expressed in one class but not the other.

## Super class

[`mlr3filters::Filter`](https://mlr3filters.mlr-org.com/reference/Filter.html)
-\> `FilterZeroProp`

## Methods

### Public methods

- [`FilterZeroProp$new()`](#method-FilterZeroProp-new)

- [`FilterZeroProp$clone()`](#method-FilterZeroProp-clone)

Inherited methods

- [`mlr3filters::Filter$calculate()`](https://mlr3filters.mlr-org.com/reference/Filter.html#method-calculate)
- [`mlr3filters::Filter$format()`](https://mlr3filters.mlr-org.com/reference/Filter.html#method-format)
- [`mlr3filters::Filter$help()`](https://mlr3filters.mlr-org.com/reference/Filter.html#method-help)
- [`mlr3filters::Filter$print()`](https://mlr3filters.mlr-org.com/reference/Filter.html#method-print)

------------------------------------------------------------------------

### Method `new()`

Creates a new instance of this Filter.

#### Usage

    FilterZeroProp$new()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    FilterZeroProp$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
