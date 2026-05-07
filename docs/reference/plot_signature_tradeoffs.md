# Plot Signature Selection Trade-offs

Visualize the trade-offs between performance, stability, and parsimony.
Shows all candidates with the selected/Pareto-optimal highlighted.

## Usage

``` r
plot_signature_tradeoffs(
  selection_result,
  x_var = "stability",
  y_var = "mean_metric",
  color_var = "mean_k"
)
```

## Arguments

- selection_result:

  Result from select_best_signature()

- x_var:

  Variable for x-axis: "mean_metric", "stability", or "mean_k"

- y_var:

  Variable for y-axis: "mean_metric", "stability", or "mean_k"

- color_var:

  Variable for color: "mean_metric", "stability", "mean_k", or
  "selected"

## Value

A ggplot2 object
