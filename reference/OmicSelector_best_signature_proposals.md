# OmicSelector_best_signature_proposals

Propose the best signture based on benchamrk methods. This function
calculated the \`metaindex\` value which is the harmonic mean of
accuracy on train, test and validation dataset. In the next step, it
sorts the miRNA sets based on \`metaIndex1\` score. The first row in
resulting data frame is the winner miRNA set.

## Usage

``` r
OmicSelector_best_signature_proposals(
  benchmark_csv = "benchmark1578929876.21765.csv",
  without_train = F
)
```

## Arguments

- benchmark_csv:

  Path to benchmark csv.

- without_train:

  One can argue, that accuracy on training dataset should not be used in
  calculation of metaIndex. By setting this to TRUE, you can calculate
  it without train dataset.

## Value

The benchmark sorted by metaIndex. First row is the best performing
miRNA set.
