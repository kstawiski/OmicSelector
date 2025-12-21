# OmicSelector_heatmap

Draw a heatmap of selected miRNAs.

## Usage

``` r
OmicSelector_heatmap(
  x = trainx[, 1:10],
  rlab = data.frame(Batch = dane$Batch, Class = dane$Class),
  zscore = F,
  margins = c(10, 10),
  expression_name = "log10(TPM)",
  trim_min = NULL,
  trim_max = NULL,
  centered_on = NULL,
  legend_pos = "topright",
  legend_cex = 0.8,
  ...
)
```

## Arguments

- x:

  Matrix of log-transformed TPM-normalized counts with miRNAs in columns
  and cases in rows.

- rlab:

  Data frame of factors to be marked on heatmap (like batch or class).
  Maximum of 2 levels for every variable is supported.

- zscore:

  Whether to z-score values before clustering and plotting.

- expression_name:

  What should be written on the plot?

- trim_min:

  Trim lower than.. Useful for setting appropriate scale

- trim_max:

  Trim greater than.. Useful for setting appropriate scale

- centered_on:

  On which value should the scale be centered? If null - median will be
  used.

- legend_pos:

  Where should the legend should be? Default: topright

- legend_cex:

  How large should the legend should be? Default: 0.8

## Value

Heatmap.
