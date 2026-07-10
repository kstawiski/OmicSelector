# Single-sample method-bank registry

Returns the package-facing implementation map for the OmicSelector
single-sample method bank. The registry is intentionally executable
metadata: each row names the exported fit, predict, score, or helper
function that implements the single-sample method or auxiliary protocol.

## Usage

``` r
singlesample_method_bank()
```

## Value

A data frame with one row per method or auxiliary protocol.
