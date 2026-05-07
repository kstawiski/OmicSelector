# Create Report Data Schema

Creates a standardized data structure containing all information
required for TRIPOD+AI compliant reporting. This separates computation
from rendering.

## Usage

``` r
create_report_data(
  benchmark_result,
  stability = NULL,
  pipeline = NULL,
  config = NULL
)
```

## Arguments

- benchmark_result:

  Result from BenchmarkService\$run()

- stability:

  Optional NogueiraStability object

- pipeline:

  OmicPipeline object used in analysis

- config:

  Optional configuration list for provenance

## Value

A ReportData object (list) with standardized structure
