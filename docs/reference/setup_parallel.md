# Configure Parallelization for OmicSelector

Sets up parallel processing using the future package. Call this before
running benchmarks or stability computations.

## Usage

``` r
setup_parallel(
  workers = NULL,
  plan = c("multisession", "multicore", "sequential"),
  memory_limit = 500 * 1024^2
)
```

## Arguments

- workers:

  Number of workers (cores) to use. Default is available cores - 1.

- plan:

  Parallelization plan: "multisession" (default, works everywhere),
  "multicore" (Unix only, more memory efficient), "sequential" (no
  parallelism)

- memory_limit:

  Per-worker memory limit in bytes (future.globals.maxSize)

## Value

Invisible previous plan (for restoration)

## Examples

``` r
if (FALSE) { # \dontrun{
# Use 4 cores
setup_parallel(workers = 4)

# Run your analysis
result <- service$run(parallel = TRUE)

# Reset to sequential
reset_parallel()
} # }
```
