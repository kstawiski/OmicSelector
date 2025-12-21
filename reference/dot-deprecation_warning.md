# Emit deprecation warning for leaky functions

Issues a prominent warning when deprecated/leaky functions are called

## Usage

``` r
.deprecation_warning(func_name, replacement = NULL, reason = NULL)
```

## Arguments

- func_name:

  Function name

- replacement:

  Suggested replacement function

- reason:

  Reason for deprecation
