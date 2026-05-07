# OmicSelector 2.4.0

**Rigorous biomarker discovery from high-dimensional omics data with zero data leakage.**

Release 2.4.0 adds the within-sample compositional scoring layer, frozen-reference denoising add-ons, qPCR non-detect imputation, the matched-null benchmark, and the provenance pre-flight gate that are benchmarked in the OmicSelector paper.

[![R-CMD-check](https://github.com/kstawiski/OmicSelector/workflows/R-CMD-check/badge.svg)](https://github.com/kstawiski/OmicSelector/actions)

---

## Validation Status

**All 11 core modules validated on TCGA pan-cancer miRNA data (10,366 samples, 2,566 features):**

| Module | Status | Key Metrics |
|--------|--------|-------------|
| OmicPipeline | PASS | Quick validation AUC: 0.936 |
| BenchmarkService | PASS | Nested CV AUC: 0.876 (honest estimate) |
| GOF Filters | PASS | Discovered miR-139/183/145 family |
| Bayesian Tuning | PASS | Training AUC: 0.984 |
| AutoXAI | PASS | 9 correlation warnings flagged |
| Stability Ensemble | PASS | 13 features at 100% stability |
| Sequential Selector (HSFS) | PASS | AUC: 0.965 |
| Synthetic Data (SMOTE) | PASS | Quality score validated |
| Calibration | PASS | Platt scaling: 0-0.56 → 0.003-0.96 |
| FrozenComBat | PASS | Batch correction validated |
| Multi-Omics | PASS | Combined AUC: 0.827 |
| Deep Learning | SKIPPED | torch not installed |

**Top Biomarkers**: miR-183-5p, miR-145-5p, miR-182-5p are established cancer biomarkers with known oncogenic/tumor suppressor roles.

---

## Overview

OmicSelector is an R package for biomarker discovery that enforces methodologically sound machine learning practices. Built on the `mlr3` ecosystem, it guarantees:

- **Zero Data Leakage**: Feature selection occurs strictly inside cross-validation folds
- **Proper Nested CV**: Separation of outer (evaluation) and inner (selection) loops
- **Feature Stability**: Nogueira Stability Index for reproducible biomarker sets
- **Multi-Objective Selection**: Balance performance, stability, and parsimony

## Paper Artifacts

The OmicSelector paper release is anchored to the package repository at
<https://github.com/kstawiski/OmicSelector>. Release 2.4.0 exposes the
paper-facing methods as exported, documented R functions:

1. Five within-sample compositional scoring methods:
   [`ws_rclr_trimmed`](man/ws_rclr_trimmed.Rd),
   [`ws_balance_ilr`](man/ws_balance_ilr.Rd),
   [`ws_alr_pivot`](man/ws_alr_pivot.Rd),
   [`ws_mad_logratio`](man/ws_mad_logratio.Rd), and
   [`ws_dominance_score`](man/ws_dominance_score.Rd).
2. Frozen-reference denoising add-ons:
   [`fit_frozen_ruv`](man/fit_frozen_ruv.Rd) /
   [`apply_frozen_ruv`](man/apply_frozen_ruv.Rd),
   [`fit_robust_pca_residual`](man/fit_robust_pca_residual.Rd) /
   [`apply_robust_pca_residual`](man/apply_robust_pca_residual.Rd), and
   [`hemolysis_index_blondal`](man/hemolysis_index_blondal.Rd) with
   [`fit_hemolysis_rr`](man/fit_hemolysis_rr.Rd) /
   [`apply_hemolysis_rr`](man/apply_hemolysis_rr.Rd).
3. qPCR non-detect handling:
   [`qpcr_nondetect_impute`](man/qpcr_nondetect_impute.Rd) wraps
   Bioconductor `nondetects`, with
   [`qpcr_nondetect_lod_fallback`](man/qpcr_nondetect_lod_fallback.Rd) as
   the documented Ct-limit fallback.
4. Matched-null benchmark:
   [`paper3_matched_null_benchmark`](man/paper3_matched_null_benchmark.Rd)
   and
   [`paper3_matched_null_benchmark_cv`](man/paper3_matched_null_benchmark_cv.Rd)
   implement stratified random-panel benchmarking with default `K = 10000`.
5. Provenance pre-flight:
   [`os_provenance_preflight`](man/os_provenance_preflight.Rd) audits a
   manifest TSV and returns `NO_OVERLAP`, `KNOWN_OVERLAP`, or
   `UNKNOWN_ACCESSION`.
6. Dependent-p combination:
   [`paper3_bh_fdr_correct_blocked`](man/paper3_bh_fdr_correct_blocked.Rd)
   reports both conservative and textbook Brown 1975 / Kost-McDermott 2002
   reference distributions.
7. Block-aware FDR:
   [`paper3_bh_fdr_correct_blocked`](man/paper3_bh_fdr_correct_blocked.Rd)
   and
   [`paper3_bh_fdr_correct_matched_null`](man/paper3_bh_fdr_correct_matched_null.Rd)
   provide the block-aware and matched-null BH correction helpers.
8. AUC confidence intervals:
   [`paper3_hanley_mcneil_auc_ci`](man/paper3_hanley_mcneil_auc_ci.Rd)
   exposes the Hanley-McNeil 1982 95% AUC CI calculation used by the
   matched-null benchmark outputs.
9. ILR coverage gate:
   [`ws_balance_ilr`](man/ws_balance_ilr.Rd) exposes
   `min_balance_coverage` so panels that preserve the fixed balance
   dictionary can activate ILR scoring at the documented coverage threshold.

## Within-Sample Biomarker Panel Classification

OmicSelector 2.4.0 adds within-sample classification methods for
reference-cohort-free biomarker panels and provenance-aware validation helpers:

- Four image encodings: `encode_simple_grid()`, `encode_corr_grid()`,
  `encode_deepinsight()`, and `encode_ratio_image()`
- Two dense CoDA learners: `train_ratio_cnn()` for pairwise log-ratio images and
  `train_clr_mlp()` for CLR-transformed vectors
- A sparse balance interface: `train_codacore()`, which uses the official
  `codacore` backend when available and otherwise falls back to a deterministic
  sparse-balance model
- `mlr3` learners registered on load: `classif.ratio_cnn`, `classif.clr_mlp`,
  and `classif.codacore`
- Group-aware resampling, provenance-floor, identifiability, matched-null, and
  operating-point gates for public-omics biomarker-panel audits

```r
library(OmicSelector)

X <- matrix(rnorm(14 * 40), nrow = 40, ncol = 14)
y <- factor(rep(c("control", "case"), each = 20), levels = c("control", "case"))

clr_fit <- train_clr_mlp(X, y, backend = "glm", verbose = FALSE)
ratio_fit <- train_ratio_cnn(
  X_train = X,
  y_train = y,
  X_test = X[1:5, , drop = FALSE],
  epochs = 5,
  verbose = FALSE
)

mlr3::lrn("classif.clr_mlp", backend = "glm", verbose = FALSE)
```

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

## Web UI (Shiny)

Run the interactive WebUI locally:

```r
shiny::runApp("shiny/omicselector", host = "0.0.0.0", port = 3838)
```

Run via Docker (Shiny WebUI + RStudio Server):

```bash
docker compose up
```

Then open:
- WebUI: `http://localhost:3838`
- RStudio Server: `http://localhost:8787` (user: `rstudio`, password: `omicselector`)

Uploads and results are saved under `/OmicSelector/analyses/<analysis_id>` inside
the container and persisted to `./analyses` on the host. Set
`OMICSELECTOR_ANALYSES_DIR` to override the storage location.

---

## Recipes

### Binary Classification (Single-Omics)

```r
pipeline <- OmicPipeline$new(
  data = my_data,
  target = "outcome",
  positive = "Case"
)

learner <- pipeline$create_graph_learner(
  filter = "anova",
  model = "rpart",
  n_features = 20
)

fit <- pipeline$fit(learner, seed = 1)
fit$selected_features
```

### Multi-Omics (Named List Input)

```r
pipeline <- OmicPipeline$new(
  data = list(rna = rna_data, mirna = mirna_data),
  target = "outcome",
  positive = "Case"
)

learner <- pipeline$create_graph_learner(
  filter = "mrmr",
  model = "rpart",
  n_features = 30
)
```

### Benchmark (Nested CV + Screening + Caching)

```r
learner <- pipeline$create_graph_learner(
  filter = "anova",
  model = "rpart",
  n_features = 20,
  screening = TRUE,
  screening_frac = 0.2
)

result <- pipeline$benchmark(
  learners = learner,
  outer_folds = 5,
  inner_folds = 3,
  seed = 42,
  cache_dir = "cache",
  parallel = TRUE,
  threads = 1
)
```

### Deep Learning (mlr3torch + Autoencoder)

```r
# Requires: torch + mlr3torch
learner <- pipeline$create_graph_learner(
  filter = "anova",
  model = "mlp",
  n_features = 50,
  autoencoder = list(
    latent_dim = 32,
    hidden_layers = c(128, 64),
    epochs = 50,
    batch_size = 64,
    early_stopping = TRUE,
    patience = 10
  )
)

result <- pipeline$benchmark(
  learners = learner,
  outer_folds = 3,
  inner_folds = 2,
  seed = 1
)
```

Transfer learning via autoencoder pretraining:

```r
ae <- autoencoder_fit(x_train, latent_dim = 32, epochs = 50)
autoencoder_save(ae, "pretrained_ae.pt")

learner <- pipeline$create_graph_learner(
  filter = "anova",
  model = "mlp",
  n_features = 50,
  autoencoder = list(
    pretrained = "pretrained_ae.pt",
    freeze_encoder = TRUE
  )
)
```

Checkpoint export/import for fine‑tuning:

```r
fit <- pipeline$fit(learner, seed = 1)
export_omicfit_checkpoint(fit, "mlp_checkpoint.pt")

finetuned <- finetune_mlr3torch_checkpoint(
  learner = learner,
  task = pipeline$get_task(),
  path = "mlp_checkpoint.pt",
  epochs = 30
)

# Or stay at the OmicFit level
finetuned_fit <- finetune_omicfit_checkpoint(
  fit = fit,
  task = pipeline$get_task(),
  path = "mlp_checkpoint.pt",
  epochs = 30
)
```

---

## Configuration Options

### Feature Selection Methods

Use the `filter` parameter in `create_graph_learner()`:

**Univariate Statistical Tests**

| Method | Code | Description | Best For |
|--------|------|-------------|----------|
| ANOVA F-test | `"anova"` | Tests mean differences between classes | Default choice, continuous features |
| Kruskal-Wallis | `"kruskal"` | Non-parametric rank-based test | Non-normal distributions |
| Chi-Squared | `"chi_squared"` | Chi-squared test | Discrete/categorical features |

**Variance-Based**

| Method | Code | Description | Best For |
|--------|------|-------------|----------|
| Variance | `"variance"` | Removes low-variance features | Pre-filtering, removing constants |

**Correlation-Based**

| Method | Code | Description | Best For |
|--------|------|-------------|----------|
| Correlation | `"correlation"` | Correlation with target | Quick univariate filter |

**Information-Theoretic**

| Method | Code | Description | Best For |
|--------|------|-------------|----------|
| Information Gain | `"information_gain"` | Entropy-based importance | Mixed feature types |
| Gain Ratio | `"gain_ratio"` | Normalized information gain | Avoiding bias toward high-cardinality |
| mRMR | `"mrmr"` | Min Redundancy Max Relevance | Reducing feature redundancy |
| CMIM | `"cmim"` | Conditional Mutual Info Max | Complex feature dependencies |
| JMIM | `"jmim"` | Joint Mutual Info Max | Capturing feature interactions |
| JMI | `"jmi"` | Joint Mutual Information | Similar to JMIM |

**Model-Based Importance**

| Method | Code | Description | Best For |
|--------|------|-------------|----------|
| AUC | `"auc"` | AUC of univariate models | Classification performance |
| Relief | `"relief"` | Instance-based algorithm | Detecting interactions |
| RF Importance | `"importance"` | Random Forest importance | Non-linear relationships |
| Permutation | `"permutation"` | Permutation importance | Model-agnostic importance |

```r
# Example: Try different filters
learner_anova <- pipeline$create_graph_learner(filter = "anova", model = "ranger", n_features = 20)
learner_mrmr <- pipeline$create_graph_learner(filter = "mrmr", model = "ranger", n_features = 20)
learner_relief <- pipeline$create_graph_learner(filter = "relief", model = "xgboost", n_features = 15)
```

### Classification Models

Use the `model` parameter in `create_graph_learner()`:

**Tree-Based Ensembles**

| Model | Code | Description | Strengths |
|-------|------|-------------|-----------|
| Random Forest | `"ranger"` | Fast RF implementation | Handles interactions, robust, default choice |
| XGBoost | `"xgboost"` | Gradient boosting | High performance, handles missing values |
| LightGBM | `"lightgbm"` | Fast gradient boosting | Very fast, memory efficient |
| Decision Tree | `"rpart"` | Single CART tree | Interpretable, fast |

**Linear Models**

| Model | Code | Description | Strengths |
|-------|------|-------------|-----------|
| Elastic Net | `"glmnet"` | L1/L2 regularized regression | Interpretable coefficients, feature selection |
| Logistic Regression | `"log_reg"` | Simple logistic regression | Baseline, highly interpretable |

**Distance-Based Methods**

| Model | Code | Description | Strengths |
|-------|------|-------------|-----------|
| SVM | `"svm"` | Support Vector Machine (RBF) | High-dimensional data |
| k-NN | `"kknn"` | K-Nearest Neighbors | Non-parametric, simple |

**Probabilistic Classifiers**

| Model | Code | Description | Strengths |
|-------|------|-------------|-----------|
| Naive Bayes | `"naive_bayes"` | Probabilistic classifier | Fast, works with small data |
| LDA | `"lda"` | Linear Discriminant Analysis | Dimensionality reduction |
| QDA | `"qda"` | Quadratic Discriminant Analysis | Non-linear decision boundaries |

**Neural Networks**

| Model | Code | Description | Strengths |
|-------|------|-------------|-----------|
| Neural Net | `"nnet"` | Single-layer neural network | Non-linear relationships |

```r
# Example: Compare models
learner_rf <- pipeline$create_graph_learner(filter = "anova", model = "ranger", n_features = 15)
learner_xgb <- pipeline$create_graph_learner(filter = "anova", model = "xgboost", n_features = 15)
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
  impute_method = "median",   # Handle missing values: "median", "mean", "sample", "pmm"
  scale = TRUE,               # Scale features to zero mean, unit variance

  # Class Imbalance (applied INSIDE each CV fold)
  oversample = NULL,          # NULL (none), "smote", or "rose"

  # Batch Effect Correction
  batch_correct = FALSE       # TRUE to apply FrozenComBat (requires batch column)
)
```

#### Advanced Imputation: PMM (Predictive Mean Matching)

For complex missing data patterns, use PMM imputation via the `mice` package:

```r
# PMM uses predictive models to impute missing values
# Best for: MAR (Missing At Random) data, preserving distributions
learner <- pipeline$create_graph_learner(
  filter = "anova",
  model = "ranger",
  n_features = 20,
  impute_method = "pmm"       # Predictive Mean Matching via mice
)

# PMM advantages:
# - Preserves data distribution (unlike mean imputation)
# - Handles complex correlation structures
# - Produces plausible values from observed data
# - Applied properly inside each CV fold to prevent leakage
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

### Phase 5: Cutting-Edge Biomarker Discovery

OmicSelector includes state-of-the-art methods validated on TCGA data:

#### GOF Filters for Sparse/Zero-Inflated Data

Traditional filters (ANOVA, variance) miss biologically important features that are "off" in one condition but expressed in another. GOF filters detect these patterns:

```r
library(OmicSelector)

# Use KS filter - captures distributional differences including zero-inflation
learner <- pipeline$create_graph_learner(
  filter = "gof_ks",      # Kolmogorov-Smirnov GOF filter
  model = "ranger",
  n_features = 30
)

# Other GOF options:
# filter = "hurdle"      # Two-part model (zero-frequency + magnitude)
# filter = "zero_prop"   # Simple zero-proportion difference

# GOF filters discovered miR-200c/miR-141 family in TCGA kidney data
# that ANOVA missed - these are established EMT regulators in cancer
```

#### Bayesian Hyperparameter Optimization

Replace grid search with intelligent Bayesian optimization:

```r
# Bayesian-tuned glmnet with omics-optimized search space
autotuner <- make_autotuner_glmnet(
  task,
  n_evals = 20,           # Budget of evaluations
  inner_folds = 3         # Inner CV for tuning
)
autotuner$train(task)

# Get optimal parameters
params <- get_optimal_params(autotuner)
# Returns: alpha, s (lambda) optimized for your data
```

#### AutoXAI: Interpretability with Correlation Warnings

Get SHAP-based interpretability with automatic warnings for correlated features:

```r
# Train a model
learner <- lrn("classif.ranger", predict_type = "prob", importance = "impurity")
learner$train(task)

# Run XAI pipeline
xai <- xai_pipeline(
  learner = learner,
  task = task,
  top_k = 20,
  cor_threshold = 0.7     # Warn if features correlated > 0.7
)

# View results
print(xai$top_features)          # Top features by permutation importance
print(xai$correlations$clusters) # Correlated feature clusters
plot_xai_importance(xai)         # Visualization

# Warnings prevent SHAP misinterpretation when features are collinear
```

#### Bootstrap Stability Ensemble

Find reproducible biomarkers across resamples:

```r
# Create stability ensemble with multiple filters
ensemble <- create_stability_ensemble(
  preset = "default",     # Uses mRMR, AUC, Information Gain
  n_bootstrap = 100,      # 100 bootstrap iterations
  n_features = 30
)

# Fit and get stable features
ensemble$fit(task, seed = 42)
stable <- ensemble$get_feature_importance(30)

# Features with >90% selection frequency are highly reproducible
# TCGA validation found 6 biomarkers with 100% stability
```

#### SMOTE for Imbalanced Data

Proper SMOTE application inside CV folds:

```r
# Balance classes with SMOTE (applied inside each fold)
task_balanced <- smote_augment(task, ratio = 1.0, k = 5)

# Validate synthetic data quality
validation <- validate_synthetic(real_data, synthetic_data)
print(validation$quality_score)  # 0-1, higher is better
```

### Feature Stability Analysis

```r
# Compute Nogueira Stability Index with Bootstrap Confidence Intervals
stability <- compute_stability_from_resample(result$benchmark_result)
print(stability)
#> Nogueira Stability Index: 0.823
#> Interpretation: Good stability

# Bootstrap CI for uncertainty quantification
stability_ci <- compute_nogueira_stability(
  feature_matrix,
  n_bootstrap = 1000,      # Number of bootstrap iterations
  conf_level = 0.95        # 95% confidence interval
)
print(stability_ci)
#> Nogueira Stability: 0.823
#> 95% Bootstrap CI: [0.756, 0.889]
#> Interpretation: Good stability (CI does not include 0.7 threshold)

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

# Apply to new predictions (calibrator is a function)
calibrated_probs <- calibrator(new_probabilities)
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

### Quality Control Validation

OmicSelector includes rigorous QC checks for clinical-grade biomarker discovery:

```r
# Run comprehensive QC validation
qc_results <- project$run_qc_validation()

# QC checks include:
# - Sample size adequacy (minimum samples per class)
# - Dimensionality ratio (features vs samples)
# - Events Per Variable (EPV) for model stability
# - Missing data patterns
# - Feature quality (variance, zero-inflation)

print(qc_results)
#> QC Validation Results:
#>   sample_size:    PASS (n=200, min_class=87)
#>   dimensionality: WARN (p/n ratio = 2.5, high dimensionality)
#>   epv:            PASS (EPV = 8.7, adequate for selected features)
#>   missing_data:   PASS (0.3% missing, handled with PMM)
#>   feature_quality: PASS (no constant features detected)

# Warnings guide corrective action without blocking workflow
```

### Model Card Generation

Generate clinical-ready Model Cards for stakeholder communication:

```r
# After training a final model
model_card <- project$generate_model_card(
  model_name = "Glioblastoma miRNA Classifier",
  intended_use = "Diagnostic aid for glioblastoma classification",
  target_population = "Adult patients with suspected CNS tumors",
  clinical_context = "Pre-surgical tumor characterization"
)

# Model Card includes:
# - Performance metrics with confidence intervals
# - Training data characteristics
# - Intended use and limitations
# - Feature stability assessment
# - TRIPOD+AI compliance notes

# Export as HTML for stakeholders
export_model_card(model_card, "model_card.html")
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

### Core Pipeline

| Module | Description |
|--------|-------------|
| `OmicPipeline` | Build mlr3 graphs with zero-leakage preprocessing |
| `BenchmarkService` | Nested CV with proper inner/outer loop separation |
| `select_best_signature()` | Multi-objective signature selection |
| `compute_nogueira_stability()` | Feature selection stability metrics |

### Batch & Calibration

| Module | Description |
|--------|-------------|
| `FrozenComBat` | Batch correction with frozen parameters |
| `fit_platt_scaling()` | Platt calibration |
| `fit_isotonic_calibration()` | Isotonic regression calibration |

### Phase 5: Advanced Features

| Module | Description |
|--------|-------------|
| `FilterGOF_KS` | Kolmogorov-Smirnov filter for sparse data |
| `FilterHurdle` | Two-part hurdle model filter |
| `FilterZeroProp` | Zero-proportion difference filter |
| `make_autotuner_glmnet()` | Bayesian-tuned elastic net |
| `make_autotuner_xgboost()` | Bayesian-tuned XGBoost |
| `make_autotuner_ranger()` | Bayesian-tuned Random Forest |
| `xai_pipeline()` | DALEX-based interpretability with correlation warnings |
| `plot_xai_importance()` | Feature importance visualization |
| `create_stability_ensemble()` | Bootstrap stability feature selection |
| `smote_augment()` | SMOTE for class imbalance |
| `validate_synthetic()` | Synthetic data quality validation |

### Quality Control & Clinical Reporting

| Module | Description |
|--------|-------------|
| `validate_data_quality()` | Sample size, dimensionality, EPV checks |
| `generate_model_card()` | Clinical-ready Model Card for stakeholders |
| `generate_tripod_report()` | TRIPOD+AI compliant reporting |

### Multi-Omics & Export

| Module | Description |
|--------|-------------|
| `MultiOmicsStacker` | Late integration of multi-omics data |
| `export_vetiver()` | Model export for deployment |
| `export_bundle()` | Self-contained model bundle with metadata |
| `export_onnx()` | ONNX export for cross-platform deployment |

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

## Publications Using OmicSelector

If you use OmicSelector in your research, please cite:

1. **Stawiski K**, Fortner RT, Elias KM, Fendler W, on behalf of the miRPOC Consortium. *External validation of a published circulating miRNA ovarian-cancer signature with OmicSelector-guided stability analysis in 13,411 samples.* (in preparation for EBioMedicine, 2026). — Used FrozenComBat, StabilityEnsemble, GOF Filters, calibration, and SHAP modules in a zero-leakage LODO-CV framework across 8 3D-Gene microarray datasets and cross-platform transfer to FirePlex immunoassay.

2. **Stawiski K**, Fortner RT, Pestarino L, et al. *Circulating miRNAs as pre-diagnostic biomarkers for ovarian cancer detection using serum from the Norwegian Janus Serum Bank.* Front Oncol. 2024;14:1389066. PMID: 38983926.

## Authors

- [Konrad Stawiski, M.D., Ph.D.](https://konsta.com.pl) (konrad.stawiski@umed.lodz.pl)
- Marcin Kaszkowiak, M.D.
- Damian Mikulski, M.D.

Supervised by: Prof. Wojciech Fendler, M.D., Ph.D.

Department of Biostatistics and Translational Medicine, Medical University of Lodz, Poland

## Links

- [Documentation](https://kstawiski.github.io/OmicSelector/)
- [Issues & Bug Reports](https://github.com/kstawiski/OmicSelector/issues)
- [Source Code](https://github.com/kstawiski/OmicSelector)
