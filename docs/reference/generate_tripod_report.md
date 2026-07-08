# Generate TRIPOD+AI Report

Renders a TRIPOD+AI compliant report from a ReportData object.

## Usage

``` r
generate_tripod_report(
  results,
  output_file = "tripod_report.html",
  format = c("html", "pdf"),
  template = NULL
)
```

## Arguments

- results:

  Benchmark results or ReportData object

- output_file:

  Output file path

- format:

  Output format: "html" (default) or "pdf"

- template:

  Custom Rmd template path (optional)

## Value

Path to generated report (invisibly)
