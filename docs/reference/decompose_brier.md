# Decompose Brier Score

Decomposes Brier score into reliability (calibration), resolution
(discrimination), and uncertainty components.

## Usage

``` r
decompose_brier(probs, labels)
```

## Arguments

- probs:

  Numeric vector of predicted probabilities

- labels:

  Binary labels (0/1)

## Value

A list containing: - brier: Total Brier score - reliability: Calibration
component (lower is better) - resolution: Discrimination component
(higher is better) - uncertainty: Inherent uncertainty (fixed for given
class balance)

## Details

Brier score decomposes as: \$\$Brier = Reliability - Resolution +
Uncertainty\$\$

\- \*\*Reliability\*\* (calibration): measures deviation from perfect
calibration - \*\*Resolution\*\*: measures how much predictions differ
across groups - \*\*Uncertainty\*\*: base rate of positive class, p(1-p)

A well-calibrated model has low reliability and high resolution.

## Examples

``` r
if (FALSE) { # \dontrun{
probs <- c(0.1, 0.4, 0.6, 0.9)
labels <- c(0, 0, 1, 1)
decompose_brier(probs, labels)
} # }
```
