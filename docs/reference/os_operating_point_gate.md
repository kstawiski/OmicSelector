# Fail-Closed Operating-Point Gate

Fail-Closed Operating-Point Gate

## Usage

``` r
os_operating_point_gate(
  operating_points,
  identifiability_status,
  min_specificity = 0.95,
  min_sensitivity = 0.5
)
```

## Arguments

- operating_points:

  Object from `os_operating_points`.

- identifiability_status:

  Terminal status from
  [`os_identifiability_gate()`](https://kstawiski.github.io/OmicSelector/reference/os_identifiability_gate.md)
  or an equivalent provenance review.

- min_specificity:

  Minimum specificity required.

- min_sensitivity:

  Minimum sensitivity required at `min_specificity`.

## Value

A one-row data frame with gate status.
