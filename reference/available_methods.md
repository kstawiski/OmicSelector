# Inspect installed method dependencies

Checks installed R packages, importable Python modules and cached
weights. Does not install or download. Imports initialize the configured
Python interpreter. Dependency availability does not establish that a
particular panel or metadata set satisfies the existing method audits.

## Usage

``` r
available_methods(methods = NULL, families = NULL, python = TRUE, weights = TRUE)
```

## Arguments

- methods:

  Optional method identifiers.

- families:

  Optional roster family letters; combined with methods as a union.

- python:

  Check imports in the active reticulate environment. FALSE leaves
  Python methods unverified.

- weights:

  Check local foundation weights without downloading. FALSE leaves
  weight requirements unverified.

## Value

Requirements plus available, missing_r, missing_python, missing_weights,
reason and install_command. Unchecked or undeclared requirements give NA
availability.

## Examples

``` r
available_methods("ws-rclr-trimmed", python = FALSE)
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
#>                                                     cost_status available
#> 1 median of completed PB1 CPU cells; see measurements attribute      TRUE
#>   missing_r missing_python missing_weights               reason
#> 1                                          dependencies present
#>                                            install_command
#> 1 install_method_dependencies(methods = "ws-rclr-trimmed")
```
