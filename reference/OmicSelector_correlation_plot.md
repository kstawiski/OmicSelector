# OmicSelector_correlation_plot

Draw a correlation scatterplot with optional y = x reference and summary
stats.

## Usage

``` r
OmicSelector_correlation_plot(
  var1,
  var2,
  labvar1,
  labvar2,
  title = NULL,
  yx = TRUE,
  metoda = "pearson",
  gdzie_legenda = "topleft"
)
```

## Arguments

- var1:

  First numeric vector.

- var2:

  Second numeric vector.

- labvar1:

  Label for x-axis.

- labvar2:

  Label for y-axis.

- title:

  Optional plot title.

- yx:

  Logical, draw y = x reference line (useful for calibration plots).

- metoda:

  Correlation method: "pearson" or "spearman".

- gdzie_legenda:

  Legend position passed to \`legend()\`.
