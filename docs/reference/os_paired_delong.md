# Paired DeLong Test for AUC Comparison

Compare two prediction-score vectors on the same labels using the paired
DeLong test (Robin et al. 2011 implementation in
[`roc.test`](https://rdrr.io/pkg/pROC/man/roc.test.html)).

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
#> $auc_a
#> [1] 0.9032
#> 
#> $auc_b
#> [1] 0.4648
#> 
#> $delta
#> [1] 0.4384
#> 
#> $z
#> [1] 6.766979
#> 
#> $p_value
#> [1] 6.574941e-12
#> 
#> $ci_delta
#> [1] 0.3114234 0.5653766
#> 
#> $n_pos
#> [1] 50
#> 
#> $n_neg
#> [1] 50
#> 
#> $alternative
#> [1] "greater"
#> 
#> $method
#> [1] "DeLong"
#> 
```
