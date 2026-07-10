# Default circulating-miRNA sequential binary partition (v1)

Returns the v1 frozen sequential binary partition (SBP) for circulating
miRNA panels used by
[`ws_balance_ilr`](https://kstawiski.github.io/OmicSelector/reference/ws_balance_ilr.md).
Each entry of the returned list is a named partition with components
`numerator` and `denominator` listing miRNA features. Several partitions
use sentinel denominators `"ALL_NON_RBC"` and `"ALL_NON_PLT"` that are
resolved at runtime as "every panel feature not in the numerator," and
one partition uses the sentinel `"TOP_K_BY_ABUNDANCE"` with a `k` field,
resolved at runtime by sorting the panel.

The partition is fixed at v1 and any future revision will be tagged
`circulating_v2`, etc., with a separate justification log under the
manuscript's sequential-binary-partition justification notes.

## Usage

``` r
ws_default_sbp(version = "circulating_v1")
```

## Arguments

- version:

  Character. Currently only `"circulating_v1"` is supported.

## Value

Named list of balance partitions.
