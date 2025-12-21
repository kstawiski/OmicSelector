# OmicSelector_OmicSelector (DEPRECATED)

\*\*DEPRECATED:\*\* This function is deprecated in OmicSelector 2.0 due
to data leakage issues. It applies SMOTE oversampling BEFORE
cross-validation splits, which causes synthetic data to leak between
folds.

The function is preserved for backward compatibility but will emit
deprecation warnings. Users should migrate to \`OmicPipeline\` for
scientifically valid results.

## Usage

``` r
OmicSelector_OmicSelector_wrapper(
  wd = getwd(),
  m = c(1:70),
  max_iterations = 10,
  code_path = system.file("extdata", "", package = "OmicSelector"),
  register_parallel = TRUE,
  clx = NULL,
  stamp = as.numeric(Sys.time()),
  prefer_no_features = 11,
  conda_path = "/home/konrad/anaconda3/bin/conda",
  debug = FALSE,
  timeout_sec = 172800,
  type = "auto"
)
```

## Arguments

- wd:

  Working directory path (default: current directory)

- m:

  Numeric vector of method IDs to run (default: 1:70)

- ...:

  Additional arguments passed to legacy function

## Value

The list of selected formulas (with deprecation warning)

## Details

\## Why is this deprecated?

The legacy \`OmicSelector_OmicSelector\` function has several critical
issues:

1\. \*\*SMOTE Leakage\*\*: Oversampling is applied to the entire
training set BEFORE cross-validation. This means synthetic samples can
appear in validation folds, leading to overly optimistic performance
estimates.

2\. \*\*Filter-then-CV\*\*: Feature selection runs on the full training
set before any CV fold isolation. Selected features may be driven by
statistical artifacts that won't generalize.

3\. \*\*Global State\*\*: Uses \`setwd()\` and file-based I/O, making
results difficult to reproduce.

\## Migration Path

Replace legacy workflows with the new mlr3-based \`OmicPipeline\` class:

“\`r \# Legacy (DEPRECATED - has leakage) \#
OmicSelector_OmicSelector(wd = ".", m = c(1, 2, 3))

\# New (zero leakage by construction) pipeline \<- OmicPipeline\$new(
data = my_data, target = "Class", positive = "Case" )

learner \<- pipeline\$create_graph_learner( filter = "anova", model =
"ranger", n_features = 20, oversample = "smote" \# Applied inside CV
folds )

service \<- BenchmarkService\$new(pipeline)
service\$add_learner(learner) result \<- service\$run() “\`

## See also

\[OmicPipeline\], \[BenchmarkService\], \[list_deprecated_functions()\]
