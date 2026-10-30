# Requirements and measured costs for registered methods

Lists dependencies beyond core Imports and, when available, measured PB1
fit time, per-specimen scoring time, peak memory and frozen model size.
Measurements describe the recorded data, hardware and settings, not a
speed guarantee.

## Usage

``` r
method_requirements(methods = NULL, families = NULL)
```

## Arguments

- methods:

  Optional method identifiers. NULL selects all unless families is set.

- families:

  Optional roster family letters. Methods and families select their
  union.

## Value

A data frame with one row per method. Dependencies are separated by
semicolons. Optional R implementations are listed separately. Cost
columns are medians of completed CPU cells. Missing measurements are NA.
The `measurements` attribute contains full cells with source and
hardware. User callbacks without dependency declarations have
`requirements_known = FALSE`.

## Examples

``` r
method_requirements("ws-rclr-trimmed")
#>         method_id family tier dep_route                           family_name
#> 1 ws-rclr-trimmed      A   R1    base-R Compositional within-sample baselines
#>   r_packages optional_r_packages python_modules model_weights      gpu
#> 1                                                             not used
#>                                       requirements_source requirements_known
#> 1 v2.6.9 registry and dependency guards; PB1 source audit               TRUE
#>   fit_seconds score_seconds_per_specimen peak_memory_bytes frozen_model_bytes
#> 1   0.3648827                 0.02434671         178630656              72429
#>                                                                               measuring_node
#> 1 Intel(R) Xeon(R) Platinum 8380 CPU @ 2.30GHz: Intel(R) Xeon(R) Platinum 8380 CPU @ 2.30GHz
#>                                                     cost_status
#> 1 median of completed PB1 CPU cells; see measurements attribute
```
