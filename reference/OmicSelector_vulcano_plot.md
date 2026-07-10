# OmicSelector_vulcano_plot

Draw a volcano plot of selected features from a differential expression
table.

## Usage

``` r
OmicSelector_vulcano_plot(
  selected_miRNAs,
  DE,
  only_label = NULL,
  take_adjusted_p = FALSE,
  log2fc_col = "log2FC",
  p_col = "p-value",
  p_adj_col = "p-value BH",
  title = NULL
)
```

## Arguments

- selected_miRNAs:

  Character vector of selected feature names.

- DE:

  Differential expression data.frame containing columns for log2FC and
  p-values.

- only_label:

  Optional character vector: label only these features (must be in
  DE\$miR).

- take_adjusted_p:

  Use adjusted p-values (BH) if TRUE.

- log2fc_col:

  Column name for log2 fold-change (default: "log2FC").

- p_col:

  Column name for raw p-values (default: "p-value").

- p_adj_col:

  Column name for adjusted p-values (default: "p-value BH").

- title:

  Optional plot title.

## Value

A ggplot object.
