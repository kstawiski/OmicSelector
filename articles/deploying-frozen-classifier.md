# Deploying a frozen single-sample classifier

Fit the reference scorer on training data, save it, and score a new
specimen. The example uses simulated abundances and needs only
OmicSelector’s core R Imports. Rendering this document also needs knitr
and rmarkdown.

``` r

library(OmicSelector)
```

    ## OmicSelector 2.6.10 - Leakage-Resistant Biomarker Discovery
    ## Use OmicPipeline$new() to start. Run ?OmicSelector for help.

``` r

set.seed(2610)
X <- matrix(rexp(41 * 30), 41, 30,
            dimnames = list(NULL, paste0("marker", 1:30)))
y <- rep(0:1, each = 20)
X[21:40, 1:5] <- X[21:40, 1:5] * 3
X_train <- X[1:40, , drop = FALSE]
x_new <- X[41, , drop = FALSE]
frozen <- deploy_singlesample(X_train, y, method = "ws-rclr-trimmed")
model_file <- tempfile(fileext = ".rds")
saveRDS(frozen, model_file)
frozen <- readRDS(model_file)
score_specimen(frozen, x_new, strict = TRUE)
```

    ## [1] 2.538847
    ## attr(,"support")
    ##   n_available n_fitted full_support
    ## 1          30       30         TRUE

``` r

stopifnot(is_singlesample_deployable(frozen, X_train[1:8, , drop = FALSE]))
unlink(model_file)
```

The score is a signed panel score, not a calibrated disease probability.
Training fixes the features and weights. Scoring uses the incoming
specimen and that frozen state. Keep the feature names, measurement
units, assay conversion and non-detect policy consistent with training.
Use the existing frozen preprocessor for counts or Ct measurements when
needed.

The audit checks that a specimen’s score does not depend on companion
specimens and that scoring leaves the model unchanged. It does not
establish clinical performance or a diagnostic threshold.

``` r

method_requirements("ws-rclr-trimmed")
```

    ##         method_id family tier dep_route                           family_name
    ## 1 ws-rclr-trimmed      A   R1    base-R Compositional within-sample baselines
    ##   r_packages optional_r_packages python_modules model_weights      gpu
    ## 1                                                             not used
    ##                                       requirements_source requirements_known
    ## 1 v2.6.9 registry and dependency guards; PB1 source audit               TRUE
    ##   fit_seconds score_seconds_per_specimen peak_memory_bytes frozen_model_bytes
    ## 1   0.3648827                 0.02434671         178630656              72429
    ##                                                                               measuring_node
    ## 1 Intel(R) Xeon(R) Platinum 8380 CPU @ 2.30GHz: Intel(R) Xeon(R) Platinum 8380 CPU @ 2.30GHz
    ##                                                     cost_status
    ## 1 median of completed PB1 CPU cells; see measurements attribute

``` r

available_methods("ws-rclr-trimmed", python = FALSE)
```

    ##         method_id family tier dep_route                           family_name
    ## 1 ws-rclr-trimmed      A   R1    base-R Compositional within-sample baselines
    ##   r_packages optional_r_packages python_modules model_weights      gpu
    ## 1                                                             not used
    ##                                       requirements_source requirements_known
    ## 1 v2.6.9 registry and dependency guards; PB1 source audit               TRUE
    ##   fit_seconds score_seconds_per_specimen peak_memory_bytes frozen_model_bytes
    ## 1   0.3648827                 0.02434671         178630656              72429
    ##                                                                               measuring_node
    ## 1 Intel(R) Xeon(R) Platinum 8380 CPU @ 2.30GHz: Intel(R) Xeon(R) Platinum 8380 CPU @ 2.30GHz
    ##                                                     cost_status available
    ## 1 median of completed PB1 CPU cells; see measurements attribute      TRUE
    ##   missing_r missing_python missing_weights               reason
    ## 1                                          dependencies present
    ##                                            install_command
    ## 1 install_method_dependencies(methods = "ws-rclr-trimmed")

``` r

install_method_dependencies(c("ai-tabpfn", "tab-tabicl"), dry_run = TRUE)
```

    ## $methods
    ## [1] "ai-tabpfn"  "tab-tabicl"
    ## 
    ## $r_packages
    ## [1] "reticulate"
    ## 
    ## $missing_r
    ## [1] "reticulate"
    ## 
    ## $python_packages
    ## [1] "numpy"           "torch"           "tabpfn==8.0.8"   "huggingface-hub"
    ## [5] "tabicl==2.1.1"  
    ## 
    ## $weights
    ## [1] "tabpfn-v2-classifier-finetuned-zk73skhh.ckpt"
    ## [2] "tabicl-classifier-v2-20260212.ckpt"          
    ## 
    ## $envname
    ## [1] "omicselector"
    ## 
    ## $lib
    ## [1] "/home/runner/.cache/R/renv/library/OmicSelector-db2f40ba/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"
    ## 
    ## $dry_run
    ## [1] TRUE

Use
[`method_requirements()`](https://kstawiski.github.io/OmicSelector/reference/method_requirements.md)
to compare required software and recorded PB1 costs. The `measurements`
attribute retains the individual datasets, measuring node, CPU, device
and settings. Missing timings mean unmeasured or unsuccessful; they do
not mean zero cost. Runtime also depends on panel and training-set size.

To add a method, run its printed `install_command`. Run the Python
installer in a fresh R session. In later sessions set
`RETICULATE_PYTHON` to the virtual environment’s interpreter before
loading OmicSelector. Python-trained methods that export plain R
parameters can retain their existing R scoring routes. Foundation models
need their recorded Python backend and weights at deployment.

The benchmark found no method with a statistically significant
improvement over the reference under the prespecified comparison. This
does not establish equal performance between methods. Choose a method
that fits the assay, available metadata and deployment resources, then
validate its performance for the intended setting. The release’s
`paper/` reproduction folder records the analysis and source tables when
the completed R16 bundle is included.
