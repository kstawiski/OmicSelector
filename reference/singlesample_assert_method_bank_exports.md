# Assert that single-sample method-bank functions are present and exported

Assert that single-sample method-bank functions are present and exported

## Usage

``` r
singlesample_assert_method_bank_exports(
  bank = singlesample_method_bank(),
  allow_deferred = TRUE
)
```

## Arguments

- bank:

  Optional registry returned by \[singlesample_method_bank()\].

- allow_deferred:

  Logical. If \`TRUE\`, deferred methods such as TabPFN-v2 are not
  treated as failures.

## Value

Invisibly returns the registry if all required functions are exported.
