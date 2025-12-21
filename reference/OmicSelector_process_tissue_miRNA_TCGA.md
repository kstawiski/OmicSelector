# OmicSelector_process_tissue_miRNA_TCGA

Process the data downloaded from TCGA.

## Usage

``` r
OmicSelector_process_tissue_miRNA_TCGA(
  data_folder = getwd(),
  remove_miRNAs_with_null_var = T
)
```

## Arguments

- data_folder:

  Directory where TCGA data were downloaded.

- remove_miRNAs_with_null_var:

  Wheter to remove the miRNAs without any expression? Default: True
