# Hanley-McNeil 1982 AUC confidence interval

Computes the closed-form Hanley-McNeil 1982 standard error and Wald
confidence interval for a binary-classification AUC from the AUC
estimate and the number of positive and negative samples. This is the
single-sample scoring bank's default per-cell AUC interval helper; it is
conservative for cross-validated estimates because it does not model
fold dependence.

## Usage

``` r
singlesample_hanley_mcneil_auc_ci(auc, n_pos, n_neg, conf_level = 0.95)
```

## Arguments

- auc:

  Numeric scalar in \[0, 1\].

- n_pos:

  Integer count of positive/case samples.

- n_neg:

  Integer count of negative/control samples.

- conf_level:

  Numeric confidence level in (0, 1). Default 0.95.

## Value

Named numeric vector with `ci_lo`, `ci_hi`, and `se`. Returns `NA`
values when inputs are not estimable.

## References

Hanley JA, McNeil BJ. (1982) The meaning and use of the area under a
receiver operating characteristic (ROC) curve. *Radiology* 143(1):
29-36.

## Examples

``` r
singlesample_hanley_mcneil_auc_ci(auc = 0.72, n_pos = 50, n_neg = 50)
#>      ci_lo      ci_hi         se 
#> 0.62007917 0.81992083 0.05098095 
```
