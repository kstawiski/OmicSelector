# OmicSelector 2.0

![](vignettes/logo.png)

**Rigorous biomarker discovery from high-dimensional omics data.**

[![R-CMD-check](https://github.com/kstawiski/OmicSelector/workflows/R-CMD-check/badge.svg)](https://github.com/kstawiski/OmicSelector/actions)

## Overview

OmicSelector is an R package for biomarker discovery that enforces methodologically sound machine learning practices. Built on the `mlr3` ecosystem, it guarantees:

- **Zero Data Leakage**: Feature selection occurs strictly inside cross-validation folds
- **Proper Nested CV**: Separation of outer (evaluation) and inner (selection) loops
- **Feature Stability**: Nogueira Stability Index for reproducible biomarker sets
- **Multi-Objective Selection**: Balance performance, stability, and parsimony

## Installation

```r
remotes::install_github("kstawiski/OmicSelector")
```

## Quick Start

```r
library(OmicSelector)
library(mlr3)

# Load your data (features + target column)
data <- read.csv("expression.csv")

# Build pipeline with embedded feature selection
pipeline <- OmicPipeline$new(data = data, target = "outcome")
graph_learner <- pipeline$create_graph_learner(
  filter = "anova",
  model = "ranger",
  n_features = 20
)

# Run nested cross-validation
benchmark <- BenchmarkService$new(
  pipeline$task,
  outer_folds = 5,
  inner_folds = 3
)
benchmark$add_learner(graph_learner)
result <- benchmark$run()

# Select optimal signature balancing performance and stability
best <- select_best_signature(result, mode = "weighted")
```

## Key Modules

| Module | Description |
|--------|-------------|
| **OmicPipeline** | Build mlr3 graphs with preprocessing |
| **BenchmarkService** | Nested CV with zero leakage |
| **select_best_signature** | Multi-objective signature selection |
| **FrozenComBat** | Batch correction with frozen parameters |
| **CalibrationService** | Probability calibration |
| **InterpretabilityService** | SHAP values with correlation warnings |
| **MultiOmicsStacker** | Late integration of multi-omics data |

## Docker

```bash
docker build -f Dockerfile.core -t omicselector:core .
docker run -it --rm -v $(pwd):/workspace omicselector:core R
```

## Citation

```
Stawiski K, Kaszkowiak M, Mikulski D, et al. OmicSelector: automatic feature
selection and deep learning modeling for omic experiments. bioRxiv. 2022.
doi: https://doi.org/10.1101/2022.06.01.494299
```

## Authors

- [Konrad Stawiski, M.D., Ph.D.](https://konsta.com.pl)
- Marcin Kaszkowiak, M.D.
- Damian Mikulski, M.D.

Supervised by: Prof. Wojciech Fendler, M.D., Ph.D.

Department of Biostatistics and Translational Medicine, Medical University of Lodz, Poland

[Issues](https://github.com/kstawiski/OmicSelector/issues) | [Documentation](https://biostat.umed.pl/OmicSelector/)
