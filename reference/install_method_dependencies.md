# Install dependencies for selected methods

Installs missing required R packages, prepares a reticulate virtual
environment and obtains cached foundation checkpoints. Repeated
installation retains existing R packages and Python packages satisfying
the plan. Foundation module versions are pinned to the inspected
environment. Other Python requirements use the selected package
resolver. The installer does not request upgrades of requirements
already satisfied.

## Usage

``` r
install_method_dependencies(methods = NULL, families = NULL, python = TRUE,
  dry_run = FALSE, envname = "omicselector", lib = .libPaths()[1L],
  repos = "https://cloud.r-project.org")
```

## Arguments

- methods:

  Optional method identifiers. NULL with NULL families selects all
  built-ins.

- families:

  Optional family letters; combined with methods as a union.

- python:

  Install Python modules and download required weights. FALSE installs
  only R requirements.

- dry_run:

  Return the plan without installation, Python imports or downloads.

- envname:

  Reticulate virtual environment name or path.

- lib:

  R library for newly installed packages.

- repos:

  CRAN repository.

## Value

A dependency plan, invisibly after installation.

## Details

When no Python interpreter is available, reticulate installs Miniconda
before creating the environment. Use a fresh R session if another Python
interpreter is initialized. Set RETICULATE_PYTHON to the returned
environment's interpreter in later sessions. The selbal package comes
from malucalle/selbal through remotes. Network access is needed for
missing packages and weights. The installer does not alter model-fitting
defaults or install optional alternative implementations.

## Examples

``` r
install_method_dependencies("ai-tabpfn", dry_run = TRUE)
#> $methods
#> [1] "ai-tabpfn"
#> 
#> $r_packages
#> [1] "reticulate"
#> 
#> $missing_r
#> [1] "reticulate"
#> 
#> $python_packages
#> [1] "numpy"           "torch"           "tabpfn==8.0.8"   "huggingface-hub"
#> 
#> $weights
#> [1] "tabpfn-v2-classifier-finetuned-zk73skhh.ckpt"
#> 
#> $envname
#> [1] "omicselector"
#> 
#> $lib
#> [1] "/home/runner/.cache/R/renv/library/OmicSelector-db2f40ba/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"
#> 
#> $dry_run
#> [1] TRUE
#> 
```
