# OmicSelector_profileplot

Profile plot for comparing score patterns across methods or samples.

## Usage

``` r
OmicSelector_profileplot(
  form,
  Method.id,
  standardize = TRUE,
  interval = 10,
  by.pattern = TRUE,
  original.names = TRUE
)
```

## Arguments

- form:

  Numeric matrix/data.frame of scores (rows = methods, cols =
  subscores).

- Method.id:

  Optional vector of method identifiers (length = nrow(form)).

- standardize:

  Logical; if TRUE, z-score each column.

- interval:

  Number of intervals on y-axis for the pattern view.

- by.pattern:

  Logical; if TRUE, draw pattern view with ggplot2.

- original.names:

  Logical; if TRUE, use column names as labels.

## Value

A ggplot object (by.pattern = TRUE) or NULL (base plot).
