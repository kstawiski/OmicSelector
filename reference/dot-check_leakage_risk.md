# Check if a function call would cause data leakage

Internal function to detect potential data leakage patterns

## Usage

``` r
.check_leakage_risk(func_name, args)
```

## Arguments

- func_name:

  Name of the function being called

- args:

  Arguments passed to the function

## Value

Logical indicating if leakage was detected
