# OmicSelector_merge_formulas

Merge and filter formulas. This function can be used to merge the
formulas\*.RDS created by different runs of
\`OmicSelector_OmicSelector()\` and filter them to keep the maximum
number of miRNAs. This may be useful in planning qPCR validation. The
result of this function is \`featureselection_formulas_final.RDS\` file,
which can be futher supplied to \`OmicSelector_benchmark()\`.

## Usage

``` r
OmicSelector_merge_formulas(wd = getwd(), max_miRNAs = 11, add = list())
```

## Arguments

- wd:

  Working directory with formulas\*.RDS files.

- max_miRNAs:

  Maximum number of miRNAs to be selected in formulas.

- add:

  List of additional sets that should be added to the formulas.

## Value

Final formulas object, also save as
\`featureselection_formulas_final.RDS\` in working directory.
