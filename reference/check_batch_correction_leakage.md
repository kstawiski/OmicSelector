# Check for Batch Correction Leakage

Diagnostic function to detect if batch correction was improperly applied
to all data together (leakage) vs. properly within CV folds.

## Usage

``` r
check_batch_correction_leakage(original_data, corrected_data, batch)
```

## Arguments

- original_data:

  Original data before any batch correction

- corrected_data:

  Data after batch correction

- batch:

  Batch vector

## Value

List with diagnostic information and warnings

## Details

This function compares variance reduction patterns. Proper frozen ComBat
applied within folds will show different patterns than leaky ComBat
applied to all data together.
