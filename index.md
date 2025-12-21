# OmicSelector 2.0

![](vignettes/logo.png)

**Rigorous biomarker discovery from high-dimensional omics data with zero data leakage.**

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

# Create pipeline from your data
pipeline <- OmicPipeline$new(
  data = my_data,          # data.frame with features + target
  target = "outcome",      # target column name
  positive = "Case"        # positive class for AUC
)

# Create learner with embedded feature selection
learner <- pipeline$create_graph_learner(
  filter = "anova",        # Feature selection: anova, mrmr, variance, correlation
  model = "ranger",        # Model: ranger, glmnet, svm, log_reg
  n_features = 20          # Number of features to select
)

# Run nested cross-validation
benchmark <- BenchmarkService$new(
  task = pipeline,
  outer_folds = 5,
  inner_folds = 3,
  seed = 42
)
benchmark$add_learner(learner)
result <- benchmark$run()

# Analyze stability and select best signature
stability <- compute_stability_from_resample(result$benchmark_result)
best <- select_best_signature(result, mode = "weighted")
```

## Data Format

Your data should be a `data.frame` with:
- **Feature columns**: Numeric values (gene expression, miRNA counts, etc.)
- **Target column**: Factor/character (classification) or numeric (regression)

```r
# Example structure:
#   gene_A  gene_B  gene_C  outcome
# 1   2.34    1.56    3.21     Case
# 2   1.12    2.89    0.45  Control
```

## Configuration Options

### Feature Selection Methods

| Method | Code | Best For |
|--------|------|----------|
| ANOVA F-test | `"anova"` | Default, continuous features |
| Kruskal-Wallis | `"kruskal"` | Non-normal distributions |
| Chi-Squared | `"chi_squared"` | Categorical features |
| Variance | `"variance"` | Pre-filtering |
| Correlation | `"correlation"` | Quick univariate |
| Information Gain | `"information_gain"` | Mixed feature types |
| Gain Ratio | `"gain_ratio"` | Avoiding cardinality bias |
| mRMR | `"mrmr"` | Reducing redundancy |
| CMIM/JMIM/JMI | `"cmim"`, `"jmim"`, `"jmi"` | Feature interactions |
| AUC | `"auc"` | Classification performance |
| Relief | `"relief"` | Detecting interactions |
| RF Importance | `"importance"` | Non-linear relationships |
| Permutation | `"permutation"` | Model-agnostic |

### Classification Models

| Model | Code | Strengths |
|-------|------|-----------|
| Random Forest | `"ranger"` | Handles interactions, robust |
| XGBoost | `"xgboost"` | High performance, handles missing values |
| LightGBM | `"lightgbm"` | Very fast, memory efficient |
| Elastic Net | `"glmnet"` | Interpretable coefficients |
| SVM | `"svm"` | High-dimensional data |
| Logistic Regression | `"log_reg"` | Baseline, interpretable |
| k-NN | `"kknn"` | Non-parametric |
| Naive Bayes | `"naive_bayes"` | Fast, small data |
| LDA/QDA | `"lda"`, `"qda"` | Dimensionality reduction |
| Neural Net | `"nnet"` | Non-linear relationships |
| Decision Tree | `"rpart"` | Interpretable |

## Key Modules

| Module | Description |
|--------|-------------|
| **OmicPipeline** | Build mlr3 graphs with preprocessing |
| **BenchmarkService** | Nested CV with zero leakage |
| **select_best_signature** | Multi-objective signature selection |
| **compute_nogueira_stability** | Feature selection stability metrics |
| **FrozenComBat** | Batch correction with frozen parameters |
| **fit_platt_scaling** | Probability calibration |
| **MultiOmicsStacker** | Late integration of multi-omics data |

## Phase 5: Advanced Features

| Module | Description |
|--------|-------------|
| **FilterGOF_KS / FilterHurdle** | GOF filters for sparse/zero-inflated data |
| **xai_pipeline** | DALEX-based interpretability with correlation warnings |
| **create_stability_ensemble** | Bootstrap stability for reproducible biomarkers |
| **make_autotuner_glmnet** | Bayesian hyperparameter optimization |
| **smote_augment** | SMOTE for class imbalance (inside CV) |

## Docker

```bash
docker build -f Dockerfile.core -t omicselector:2.0 .
docker run -it --rm -v $(pwd):/workspace omicselector:2.0 R
```

## Citation

```bibtex
@article{stawiski2022omicselector,
  title={OmicSelector: automatic feature selection and deep learning
         modeling for omic experiments},
  author={Stawiski, Konrad and Kaszkowiak, Marcin and Mikulski, Damian and others},
  journal={bioRxiv},
  year={2022},
  doi={10.1101/2022.06.01.494299}
}
```

## Authors

- [Konrad Stawiski, M.D., Ph.D.](https://konsta.com.pl) (konrad.stawiski@umed.lodz.pl)
- Marcin Kaszkowiak, M.D.
- Damian Mikulski, M.D.

Supervised by: Prof. Wojciech Fendler, M.D., Ph.D.

Department of Biostatistics and Translational Medicine, Medical University of Lodz, Poland

## Links

- [Documentation](https://biostat.umed.pl/OmicSelector/)
- [Issues & Bug Reports](https://github.com/kstawiski/OmicSelector/issues)
- [Source Code](https://github.com/kstawiski/OmicSelector)
