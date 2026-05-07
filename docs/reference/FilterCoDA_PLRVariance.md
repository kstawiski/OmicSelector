# Pairwise Log-Ratio Variance Filter

\`mlr3filters::Filter\` wrapper around \[codaFS_plr_variance()\].

## Dictionary

This Filter can be instantiated via the dictionary \`mlr_filters\` or
with the associated sugar function \`flt()\`: “\`
flt("coda_plr_variance") mlr_filters\$get("coda_plr_variance") “\`

## Super class

[`mlr3filters::Filter`](https://mlr3filters.mlr-org.com/reference/Filter.html)
-\> `FilterCoDA_PLRVariance`

## Methods

### Public methods

- [`FilterCoDA_PLRVariance$new()`](#method-FilterCoDA_PLRVariance-new)

- [`FilterCoDA_PLRVariance$clone()`](#method-FilterCoDA_PLRVariance-clone)

Inherited methods

- [`mlr3filters::Filter$calculate()`](https://mlr3filters.mlr-org.com/reference/Filter.html#method-calculate)
- [`mlr3filters::Filter$format()`](https://mlr3filters.mlr-org.com/reference/Filter.html#method-format)
- [`mlr3filters::Filter$help()`](https://mlr3filters.mlr-org.com/reference/Filter.html#method-help)
- [`mlr3filters::Filter$print()`](https://mlr3filters.mlr-org.com/reference/Filter.html#method-print)

------------------------------------------------------------------------

### Method `new()`

Creates a new instance of this Filter.

#### Usage

    FilterCoDA_PLRVariance$new()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    FilterCoDA_PLRVariance$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
