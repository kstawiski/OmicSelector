# OmicSelector_correlation_plot

Draw a correlation plot with reference f(x)=x axis and correlation
metrics.

## Usage

``` r
OmicSelector_correlation_plot(
  var1,
  var2,
  labvar1,
  labvar2,
  title,
  yx = T,
  metoda = "pearson",
  gdzie_legenda = "topleft"
)
```

## Arguments

- var1:

  First numeric variable.

- var2:

  Second numeric variable.

- labvar1:

  Label to be written on plot for fist variable.

- labvar2:

  Label to be written on plot for second variable.

- yx:

  Logical varable. TRUE if you want f(x)=x axis for reference (usefull
  for calibration plots).

- gdzie_legenda:

  Where the legend should be placed?
