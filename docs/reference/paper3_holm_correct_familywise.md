# Holm correction for family-wise method contrasts

Applies Holm step-down correction across a family-wise contrast set
(method \\\times\\ claim-state contrasts; hemolysis-injection deltas).
Per plan v0.3.3 section 5: Holm is reserved for these family-wise
comparisons, NOT for discovery-screening matched-null tests (which use
BH-FDR).

## Usage

``` r
paper3_holm_correct_familywise(p_values)
```

## Arguments

- p_values:

  Numeric vector of contrast p-values.

## Value

Numeric vector of Holm-corrected p-values.
