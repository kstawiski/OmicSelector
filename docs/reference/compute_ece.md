# Compute Expected Calibration Error (ECE)

Calculates the Expected Calibration Error, which measures how well
predicted probabilities match observed frequencies. Lower is better.

## Usage

``` r
compute_ece(probs, labels, n_bins = 10, weighting = c("samples", "uniform"))
```

## Arguments

- probs:

  Numeric vector of predicted probabilities (0-1)

- labels:

  Binary labels (0/1 or logical)

- n_bins:

  Number of bins for grouping probabilities (default: 10)

- weighting:

  How to weight bins: "samples" (default) or "uniform"

## Value

A list containing: - ece: Expected Calibration Error (0-1, lower is
better) - mce: Maximum Calibration Error - bin_data: Data frame with
per-bin statistics

## Details

ECE is computed as: \$\$ECE = \sum\_{b=1}^{B} \frac{n_b}{N}
\|accuracy_b - confidence_b\|\$\$

where B is the number of bins, n_b is samples in bin b, N is total
samples, accuracy_b is fraction of positives in bin, and confidence_b is
mean predicted probability in bin.

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulated predictions
probs <- runif(100)
labels <- rbinom(100, 1, probs)  # Well-calibrated by construction

result <- compute_ece(probs, labels)
print(result$ece)  # Should be low
} # }
```
