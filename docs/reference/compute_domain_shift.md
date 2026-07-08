# Quantify Domain Shift Between Source and Target Platforms

Computes a per-feature domain-shift summary using the Kolmogorov-Smirnov
test and Cohen's d effect size.

## Usage

``` r
compute_domain_shift(source_data, target_data)
```

## Arguments

- source_data:

  Numeric matrix with source samples in rows and features in columns.

- target_data:

  Numeric matrix with target samples in rows and features in columns.

## Value

A \`data.table\` with columns \`feature\`, \`ks_stat\`, \`ks_p\`, and
\`effect_size_d\`.

## Examples

``` r
set.seed(42)
src <- matrix(rnorm(100 * 5), 100, 5)
tgt <- matrix(rnorm(50 * 5, mean = 2), 50, 5)

compute_domain_shift(src, tgt)
#>      feature ks_stat         ks_p effect_size_d
#>       <char>   <num>        <num>         <num>
#> 1: feature_5    0.78 4.853397e-18     -2.039276
#> 2: feature_1    0.74 2.795028e-16     -1.982399
#> 3: feature_3    0.72 1.958000e-15     -1.964085
#> 4: feature_2    0.70 1.300400e-14     -2.076781
#> 5: feature_4    0.68 8.188018e-14     -2.021904
```
