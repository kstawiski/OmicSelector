# Paired DeLong Test for AUC Comparison

Compare two prediction-score vectors on the same labels using the paired
DeLong test (Robin et al. 2011 implementation in `roc.test`).

## Usage

``` r
os_paired_delong(
  y_true,
  scores_a,
  scores_b,
  alternative = c("two.sided", "greater", "less")
)
```

## Arguments

- y_true:

  Integer or logical vector of true labels (0/1).

- scores_a:

  Numeric vector of scores from pipeline A; same length as `y_true`.

- scores_b:

  Numeric vector of scores from pipeline B; same length as `y_true`.

- alternative:

  One of `"two.sided"` (default), `"greater"` (A's AUC \> B's), or
  `"less"`.

## Value

A named list with elements `auc_a`, `auc_b`, `delta` (= auc_a - auc_b),
`z`, `p_value`, `ci_delta` (95 the AUC difference; `NA` if pROC does not
return one), `n_pos`, `n_neg`, `alternative`, and `method` (= "DeLong").

## Examples

``` r
set.seed(7)
y <- rep(c(0, 1), each = 50)
a <- y + rnorm(100, sd = 0.5)            # informative
b <- rnorm(100)                          # random
os_paired_delong(y, a, b, alternative = "greater")
#> Error: os_paired_delong() requires the 'pROC' package. Install it first.
```
