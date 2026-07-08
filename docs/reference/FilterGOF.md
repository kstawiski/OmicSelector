# Goodness-of-Fit Filters for Sparse Omics Data

Implements feature selection filters designed for sparse/zero-inflated
omics data. Standard variance-based filters often discard biologically
relevant features that are "turned off" (zero) in one condition but
expressed in another.

## Details

Two filters are provided: - \*\*FilterGOF_KS\*\*: Kolmogorov-Smirnov
test comparing distributions between classes - \*\*FilterHurdle\*\*:
Two-part hurdle model testing both zero-frequency and magnitude

These filters are designed for binary classification tasks in biomarker
discovery.
