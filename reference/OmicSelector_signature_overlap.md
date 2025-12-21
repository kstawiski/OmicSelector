# OmicSelector_signature_overlap

A function to generate venn diagram and check the overlap between
formulas.

## Usage

``` r
OmicSelector_signature_overlap(
  which_formulas = c("sig", "cfs"),
  benchmark_csv = "benchmark.csv"
)
```

## Arguments

- which_formulas:

  Which formulas to check?

- benchmark_csv:

  Which benchmark to use?

## Value

Object of \`venn()\` function which can be used for plotting venn
diagram and check the overlap.
