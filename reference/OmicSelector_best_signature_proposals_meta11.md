# OmicSelector_best_signature_proposals_meta11

Propose the best signture based on benchamrk methods. This function
calculated the \`metaIndex11\` value which is the Youden-like score on
validation set (the only one that was never used in any section of the
pipeline). Formula: \`metaIndex11 = validation sensitivitiy + validation
specificity - 1\` In the next step, it sorts the miRNA sets based on
\`metaIndex11\` score. The first row in resulting data frame is the
winner miRNA set.

## Usage

``` r
OmicSelector_best_signature_proposals_meta11(
  benchmark_csv = "benchmark1578929876.21765.csv"
)
```

## Arguments

- benchmark_csv:

  Path to benchmark csv.

## Value

The benchmark sorted by metaIndex. First row is the best performing
miRNA set.
