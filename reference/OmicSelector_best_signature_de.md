# OmicSelector_best_signature_de

As a part of checkpoint, you may want to check the differential
expression of selected miRNAs. This function uses
\`OmicSelector_differential_expression_ttest()\` to check the miRNAs on
training dataset.

## Usage

``` r
OmicSelector_best_signature_de(selected_miRNAs, use_mix = F)
```

## Arguments

- selected_miRNAs:

  Vector of selected miRNAs to be checked.

- use_mix:

  By default (i.e. FALSE) we check the differential expression only on
  training dataset. If you want to check it on whole dataset (training,
  testing and validation dataset combined) set it to TRUE.

## Value

Results of differential expression.
