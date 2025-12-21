# Create Reliability Diagram Data

Generates data for plotting a reliability (calibration) diagram.

## Usage

``` r
reliability_diagram_data(probs, labels, n_bins = 10)
```

## Arguments

- probs:

  Predicted probabilities

- labels:

  True labels

- n_bins:

  Number of bins

## Value

A data frame suitable for ggplot2

## Examples

``` r
if (FALSE) { # \dontrun{
library(ggplot2)

diagram_data <- reliability_diagram_data(probs, labels)

ggplot(diagram_data, aes(x = mean_predicted, y = fraction_positive)) +
  geom_point(aes(size = n_samples)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  xlim(0, 1) + ylim(0, 1) +
  labs(x = "Mean Predicted Probability", y = "Fraction of Positives")
} # }
```
