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
# Install from GitHub
remotes::install_github("kstawiski/OmicSelector")

# Load the package
library(OmicSelector)
```

---

## Data Format

OmicSelector expects a **data.frame** with:
- **Feature columns**: Numeric values (gene expression, miRNA counts, etc.)
- **Target column**: Factor or character for classification, numeric for regression
- **Optional metadata**: patient_id, batch, etc.

### Example Data Structure

```r
# Your data should look like this:
head(my_data)
#>   gene_A  gene_B  gene_C  gene_D  gene_E  outcome
#> 1   2.34    1.56    3.21    0.89    1.45     Case
#> 2   1.12    2.89    0.45    2.11    0.78  Control
#> 3   3.45    0.23    2.67    1.34    2.90     Case
#> ...

# Required structure:
# - All feature columns must be numeric
# - Target column can be factor, character (classification) or numeric (regression)
# - No row names required (use a patient_id column if needed)
```

### Creating Test Data

```r
# Simulate expression data for 100 samples, 50 genes
set.seed(42)
n_samples <- 100
n_genes <- 50

# Create feature matrix
features <- matrix(rnorm(n_samples * n_genes), nrow = n_samples)
colnames(features) <- paste0("gene_", 1:n_genes)

# Add some signal to first 5 genes
outcome <- factor(rep(c("Case", "Control"), each = n_samples/2))
features[outcome == "Case", 1:5] <- features[outcome == "Case", 1:5] + 1.5

# Combine into data.frame
my_data <- data.frame(features, outcome = outcome)
```

### Multi-Omics Data

For multi-omics, pass a **named list** of data.frames (same samples, same row order):

```r
# RNA-seq data (100 samples x 1000 genes)
rna_data <- data.frame(matrix(rnorm(100 * 1000), nrow = 100))
names(rna_data) <- paste0("gene_", 1:1000)

# miRNA data (100 samples x 200 miRNAs)
mirna_data <- data.frame(matrix(rnorm(100 * 200), nrow = 100))
names(mirna_data) <- paste0("miR_", 1:200)

# Target vector
outcome <- factor(rep(c("Case", "Control"), each = 50))

# Create multi-omics pipeline
pipeline <- OmicPipeline$new(
  data = list(rna = rna_data, mirna = mirna_data),
  target = outcome,
  positive = "Case"
)
# Features are automatically namespaced: rna::gene_1, mirna::miR_1, etc.
```

---

## Quick Start

### Step 1: Create Pipeline

```r
library(OmicSelector)

# Create pipeline from data
pipeline <- OmicPipeline$new(
  data = my_data,          # Your data.frame
  target = "outcome",      # Name of the target column
  positive = "Case"        # Positive class for AUC calculation
)
#> OmicPipeline created: 100 samples, 50 features, target='outcome'
```

### Step 2: Create GraphLearner

```r
# Create a learner with embedded feature selection
learner <- pipeline$create_graph_learner(
  filter = "anova",        # Feature selection method
  model = "ranger",        # Classification model
  n_features = 10          # Select top 10 features
)
```

### Step 3: Run Nested Cross-Validation

```r
# Create benchmark service
benchmark <- BenchmarkService$new(
  task = pipeline,         # Your pipeline
  outer_folds = 5,         # Outer CV folds (evaluation)
  inner_folds = 3,         # Inner CV folds (selection)
  seed = 42                # For reproducibility
)

# Add learner and run
benchmark$add_learner(learner)
result <- benchmark$run()
#> Running nested cross-validation...
#> Completed in 12.3 seconds

# View results
print(result)
```

---

## Configuration Options

### Feature Selection Methods

Use the `filter` parameter in `create_graph_learner()`:

| Method | Code | Description | Best For |
|--------|------|-------------|----------|
| ANOVA F-test | `"anova"` | Tests mean differences between classes | Default choice, continuous features |
| Variance | `"variance"` | Removes low-variance features | Pre-filtering, noisy data |
| Correlation | `"correlation"` | Correlation with target | Quick univariate filter |
| Information Gain | `"information_gain"` | Entropy-based importance | Mixed feature types |
| mRMR | `"mrmr"` | Minimum Redundancy Maximum Relevance | Reducing feature redundancy |

```r
# Example: Try different filters
learner_anova <- pipeline$create_graph_learner(filter = "anova", model = "ranger", n_features = 20)
learner_mrmr <- pipeline$create_graph_learner(filter = "mrmr", model = "ranger", n_features = 20)
```

### Classification Models

Use the `model` parameter in `create_graph_learner()`:

| Model | Code | Description | Strengths |
|-------|------|-------------|-----------|
| Random Forest | `"ranger"` | Fast RF implementation | Handles interactions, robust |
| Elastic Net | `"glmnet"` | L1/L2 regularized regression | Interpretable coefficients |
| SVM | `"svm"` | Support Vector Machine | High-dimensional data |
| Logistic Regression | `"log_reg"` | Simple logistic regression | Baseline, interpretable |

```r
# Example: Compare models
learner_rf <- pipeline$create_graph_learner(filter = "anova", model = "ranger", n_features = 15)
learner_svm <- pipeline$create_graph_learner(filter = "anova", model = "svm", n_features = 15)
learner_glmnet <- pipeline$create_graph_learner(filter = "anova", model = "glmnet", n_features = 15)
```

### Using Any mlr3 Learner

OmicSelector is built on the `mlr3` ecosystem. You can use **any** mlr3 learner directly with `BenchmarkService`:

```r
library(mlr3learners)

# Get the task from the pipeline
task <- pipeline$get_task()

# Create any mlr3 learner with probability predictions
xgboost_learner <- lrn("classif.xgboost", predict_type = "prob")
nnet_learner <- lrn("classif.nnet", predict_type = "prob")

# Wrap in a graph with feature selection (optional)
library(mlr3pipelines)
library(mlr3filters)

xgb_graph <- po("filter", filter = flt("anova"), filter.nfeat = 20) %>>%
  po("learner", xgboost_learner)

# Convert to GraphLearner
xgb_graph_learner <- as_learner(xgb_graph)
xgb_graph_learner$id <- "xgboost_anova_20"

# Use in benchmark
benchmark <- BenchmarkService$new(task = pipeline, outer_folds = 5, seed = 42)
benchmark$add_learner(xgb_graph_learner)
result <- benchmark$run()
```

See the [mlr3learners documentation](https://mlr3learners.mlr-org.com/) for all available learners.

### Full Parameter Reference

```r
learner <- pipeline$create_graph_learner(
  # Feature Selection
  filter = "anova",           # Filter method (see table above)
  n_features = 20,            # Number of features to select
                              # Use decimal < 1 for proportion (e.g., 0.1 = top 10%)

  # Preprocessing
  impute_method = "median",   # Handle missing values: "median", "mean", "sample"
  scale = TRUE,               # Scale features to zero mean, unit variance

  # Class Imbalance (applied INSIDE each CV fold)
  oversample = NULL,          # NULL (none), "smote", or "rose"

  # Batch Effect Correction
  batch_correct = FALSE       # TRUE to apply FrozenComBat (requires batch column)
)
```

---

## Comparing Multiple Configurations

```r
# Create benchmark service
benchmark <- BenchmarkService$new(pipeline, outer_folds = 5, inner_folds = 3, seed = 42)

# Add multiple learners with different configurations
benchmark$add_learner(
  pipeline$create_graph_learner(filter = "anova", model = "ranger", n_features = 10),
  id = "anova_rf_10"
)
benchmark$add_learner(
  pipeline$create_graph_learner(filter = "anova", model = "ranger", n_features = 20),
  id = "anova_rf_20"
)
benchmark$add_learner(
  pipeline$create_graph_learner(filter = "mrmr", model = "glmnet", n_features = 15),
  id = "mrmr_glmnet_15"
)
benchmark$add_learner(
  pipeline$create_graph_learner(filter = "variance", model = "svm", n_features = 25),
  id = "var_svm_25"
)

# Run all configurations
result <- benchmark$run()

# Compare performance
print(result$performance)
#>           learner_id classif.auc classif.acc classif.bbrier
#> 1       anova_rf_10       0.892       0.840          0.112
#> 2       anova_rf_20       0.901       0.850          0.098
#> 3   mrmr_glmnet_15        0.885       0.830          0.121
#> 4       var_svm_25        0.878       0.820          0.134
```

---

## Advanced Features

### Feature Stability Analysis

```r
# Compute Nogueira Stability Index
stability <- compute_stability_from_resample(result$benchmark_result)
print(stability)
#> Nogueira Stability Index: 0.823
#> Interpretation: Good stability

# Get features selected in ALL folds (consensus)
consensus_features <- stability$consensus_features
print(consensus_features)
#> [1] "gene_1" "gene_3" "gene_5"

# Get selection frequency per feature
head(stability$selection_frequency)
#>   gene_1   gene_3   gene_5   gene_2   gene_4
#>     1.00     1.00     1.00     0.80     0.60
```

### Optimal Signature Selection

```r
# Select best signature balancing multiple objectives
best <- select_best_signature(
  result,
  mode = "weighted",         # Selection mode
  stability_weight = 0.3,    # Weight for stability (vs performance)
  parsimony_weight = 0.1     # Weight for fewer features
)

# Available modes:
# - "best_auc": Highest AUC only
# - "best_stability": Highest Nogueira index only
# - "weighted": Weighted combination (default)
# - "pareto": Pareto-optimal with distance-to-utopia
```

### Model Calibration

```r
# After getting predictions from your model
probabilities <- predict(model, newdata = test_data, type = "prob")
true_labels <- test_data$outcome

# Compute calibration metrics
ece <- compute_ece(probabilities, true_labels)
brier <- decompose_brier(probabilities, true_labels)

# Fit calibration (Platt scaling or isotonic regression)
calibrator <- fit_platt_scaling(probabilities, true_labels)
# or: calibrator <- fit_isotonic_calibration(probabilities, true_labels)

# Apply to new predictions
calibrated_probs <- calibrator$calibrate(new_probabilities)
```

### Batch Correction (FrozenComBat)

```r
# Include batch information in pipeline
pipeline <- OmicPipeline$new(
  data = my_data,
  target = "outcome",
  batch = "batch_id"         # Column containing batch labels
)

# Enable batch correction (parameters fit on training data only)
learner <- pipeline$create_graph_learner(
  filter = "anova",
  model = "ranger",
  n_features = 20,
  batch_correct = TRUE       # Applies FrozenComBat inside CV
)
```

### SHAP Interpretability

```r
# Create explainer for a trained model
explainer <- create_explainer(trained_learner, task)

# Compute SHAP values
shap <- shap_values(explainer, new_data)

# Check for problematic correlations
# (SHAP can be misleading when features are highly correlated)
correlations <- check_feature_correlations(
  data = my_data,
  features = selected_features,
  threshold = 0.7            # Warn if correlation > 0.7
)
```

### Multi-Omics Late Integration

```r
# Create separate pipelines for each modality
rna_pipeline <- OmicPipeline$new(rna_data, target = outcome)
mirna_pipeline <- OmicPipeline$new(mirna_data, target = outcome)

# Create stacker
stacker <- MultiOmicsStacker$new(
  modalities = list(
    rna = rna_pipeline$get_task(),
    mirna = mirna_pipeline$get_task()
  ),
  meta_learner = mlr3::lrn("classif.glmnet")
)

# Train and predict
stacker$train()
predictions <- stacker$predict(new_data)
```

---

## YAML Configuration

For reproducible batch processing:

```yaml
# config.yaml
data:
  file: "data/expression.csv"
  target: "outcome"
  positive: "Case"

pipeline:
  filter: "anova"
  model: "ranger"
  n_features: 20
  scale: true
  impute: "median"

benchmark:
  outer_folds: 5
  inner_folds: 3
  seed: 42
  stratify: true

output:
  dir: "./results"
  report_format: "html"
```

Run via CLI:
```bash
Rscript inst/bin/omicselector run --config=config.yaml
```

---

## Docker

```bash
# Build image
docker build -f Dockerfile.core -t omicselector:2.0 .

# Run interactive R session
docker run -it --rm -v $(pwd):/workspace omicselector:2.0 R

# Run analysis with config
docker run -it --rm -v $(pwd):/workspace omicselector:2.0 \
  Rscript inst/bin/omicselector run --config=/workspace/config.yaml
```

---

## Module Reference

| Module | Description |
|--------|-------------|
| `OmicPipeline` | Build mlr3 graphs with zero-leakage preprocessing |
| `BenchmarkService` | Nested CV with proper inner/outer loop separation |
| `select_best_signature()` | Multi-objective signature selection |
| `compute_nogueira_stability()` | Feature selection stability metrics |
| `FrozenComBat` | Batch correction with frozen parameters |
| `fit_platt_scaling()` | Platt calibration |
| `fit_isotonic_calibration()` | Isotonic regression calibration |
| `shap_values()` | SHAP-based interpretability |
| `MultiOmicsStacker` | Late integration of multi-omics data |
| `export_vetiver()` | Model export for deployment |

---

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
