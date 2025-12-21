# Kolmogorov-Smirnov GOF Filter

Uses the KS test D-statistic as feature importance score. Captures
differences in distribution shape, location, and zero-frequency between
two classes without assuming a parametric distribution.

## Details

The KS test compares Empirical Cumulative Distribution Functions (ECDF).
If one class has 80 a massive difference captured by the D statistic
(bounded \[0, 1\]).

Higher D values indicate better class separation.

## Dictionary

This Filter can be instantiated via the dictionary \`mlr_filters\` or
with the associated sugar function \`flt()\`: “\` flt("gof_ks")
mlr_filters\$get("gof_ks") “\`

## References

Kolmogorov, A. N. (1933). Sulla determinazione empirica di una legge di
distribuzione.

## Super class

[`mlr3filters::Filter`](https://mlr3filters.mlr-org.com/reference/Filter.html)
-\> `FilterGOF_KS`

## Methods

### Public methods

- [`FilterGOF_KS$new()`](#method-FilterGOF_KS-new)

- [`FilterGOF_KS$clone()`](#method-FilterGOF_KS-clone)

Inherited methods

- [`mlr3filters::Filter$calculate()`](https://mlr3filters.mlr-org.com/reference/Filter.html#method-calculate)
- [`mlr3filters::Filter$format()`](https://mlr3filters.mlr-org.com/reference/Filter.html#method-format)
- [`mlr3filters::Filter$help()`](https://mlr3filters.mlr-org.com/reference/Filter.html#method-help)
- [`mlr3filters::Filter$print()`](https://mlr3filters.mlr-org.com/reference/Filter.html#method-print)

------------------------------------------------------------------------

### Method `new()`

Creates a new instance of this Filter.

#### Usage

    FilterGOF_KS$new()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    FilterGOF_KS$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
