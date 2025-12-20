![](vignettes/logo.png)

# OmicSelector 2.0

[![R-CMD-check](https://github.com/kstawiski/OmicSelector/workflows/R-CMD-check/badge.svg)](https://github.com/kstawiski/OmicSelector/actions)

**OmicSelector** is an R package for rigorous biomarker discovery from high-dimensional omics data. Built on the `mlr3` ecosystem, it enforces methodologically sound machine learning practices including proper nested cross-validation, feature stability analysis, and zero data leakage.

## Key Features

- **Zero Leakage Architecture**: Feature selection occurs strictly inside cross-validation folds via `mlr3pipelines`, preventing inflated performance estimates
- **Nested Cross-Validation**: Proper inner/outer loop separation for unbiased model evaluation
- **Feature Stability Analysis**: Nogueira Stability Index to quantify reproducibility of selected biomarkers
- **Multi-Objective Signature Selection**: Choose optimal signatures balancing performance, stability, and parsimony
- **Batch Effect Correction**: FrozenComBat for training-only parameter fitting
- **Model Calibration**: Brier score, ECE, and isotonic/Platt recalibration
- **Interpretability**: SHAP values with correlation warnings via DALEX/iml integration
- **Multi-Omics Integration**: Late integration (stacking) of multiple data modalities
- **Reproducibility**: `renv` lockfiles and Docker images for deterministic builds

## Installation

```r
# Install from GitHub
remotes::install_github("kstawiski/OmicSelector")

# Load package
library(OmicSelector)
```

## Quick Start

```r
library(OmicSelector)
library(mlr3)

# Create task from expression data
task <- as_task_classif(my_data, target = "Class")

# Build pipeline with embedded feature selection (zero leakage)
pipeline <- OmicPipeline$new(task)
pipeline$add_filter("variance", cutoff = 0.8)
pipeline$add_learner("classif.ranger")
graph_learner <- pipeline$create_graph_learner()

# Run nested cross-validation
benchmark <- BenchmarkService$new(task, outer_folds = 5, inner_folds = 3)
benchmark$add_learner(graph_learner)
result <- benchmark$run()

# Select best signature (multi-objective)
best <- select_best_signature(result, mode = "weighted")
print(best)
```

## Core Components

### OmicPipeline
R6 class for building `mlr3pipelines` graphs with embedded preprocessing:
- Variance filtering, correlation filtering
- Normalization and scaling
- SMOTE/ROSE (inside CV folds)
- Graph learner creation

### BenchmarkService
Enforces proper nested cross-validation:
- Outer loop: Performance evaluation
- Inner loop: Feature selection and hyperparameter tuning
- Stores models for post-hoc analysis

### Signature Selection
Multi-objective selection with three modes:
- `constrained_1se`: 1 standard error rule with tie-breaking by stability/parsimony
- `weighted`: Weighted combination of performance, stability, parsimony
- `pareto`: Pareto frontier with distance-to-utopia selection

### Phase 4 Modules
Advanced features for PhD-level rigor:
- **FrozenComBat**: Batch correction with frozen parameters
- **Calibration**: Model probability calibration and diagnostics
- **Interpretability**: SHAP values with correlation warnings
- **Multi-Omics**: Late integration via stacking
- **Model Export**: vetiver-based deployment support

## Docker

Reproducible environment based on `rocker/r-ver:4.4.0`:

```bash
# Build core image
docker build -f Dockerfile.core -t omicselector:core .

# Run R session
docker run -it --rm -v $(pwd):/workspace omicselector:core R
```

## Citation

```
Stawiski K, Kaszkowiak M, Mikulski D, Hogendorf P, Durczynski A, Strzelczyk J, et al.
OmicSelector: automatic feature selection and deep learning modeling for omic experiments.
bioRxiv. 2022. doi: https://doi.org/10.1101/2022.06.01.494299
```

## Authors

- [Konrad Stawiski, M.D., Ph.D.](https://konsta.com.pl) (konrad.stawiski@umed.lodz.pl)
- Marcin Kaszkowiak, M.D.
- Damian Mikulski, M.D.

Supervised by: Prof. Wojciech Fendler, M.D., Ph.D.

Department of Biostatistics and Translational Medicine, Medical University of Lodz, Poland
https://biostat.umed.pl

## Issues

Report bugs and feature requests at https://github.com/kstawiski/OmicSelector/issues
