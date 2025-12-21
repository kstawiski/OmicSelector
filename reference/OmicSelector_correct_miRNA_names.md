# OmicSelector_correct_miRNA_names

Sometinmes, when using the dataset mapped to previous versions of
miRbase we may get false mismatches due to changes in terminology. This
function uses latest version of miRbase to correct all old miRNA names
to new one.

## Usage

``` r
OmicSelector_correct_miRNA_names(temp, species = "hsa", correct_dots = T)
```

## Arguments

- temp:

  Dataset with miRNA names in columns.

- correct_dots:

  Boolean variable to correct the names after correction of dots to
  hyphens. This tries to compensate the effect of make.names() function.

## Value

Corrected dataset.
