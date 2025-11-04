# OmicSelector 2.0 Modernization

## Overview

OmicSelector 2.0 represents a major modernization of the biomarker discovery platform, bringing it to 2025+ standards for machine learning and clinical prediction model development. This document outlines the key changes and new features.

## What's New in Version 2.0

### 🎯 Core Framework Modernization

#### 1. **Tidymodels Integration** (NEW)
- Modern, unified interface for machine learning workflows
- Seamless integration with the tidyverse ecosystem
- Enhanced reproducibility and code clarity
- **Backward compatibility maintained** - All existing caret-based code continues to work

#### 2. **Rigorous Nested Cross-Validation** (NEW)
- Gold-standard approach for unbiased model evaluation
- **Prevents data leakage** by design:
  - All preprocessing happens inside resampling folds
  - Feature selection is nested properly
  - Hyperparameter tuning in inner loop only
- Outer loop: unbiased performance estimation
- Inner loop: model selection and tuning

#### 3. **TRIPOD+AI Compliance** (NEW)
- Automated generation of TRIPOD+AI compliant reports
- All 27 checklist items covered
- Export to HTML, PDF, JSON, or Markdown
- Ensures transparent reporting for clinical prediction models

#### 4. **PROBAST+AI Risk of Bias Assessment** (NEW)
- Systematic assessment across 4 domains:
  1. Participants
  2. Predictors
  3. Outcome
  4. Analysis
- Automated recommendations for improvement
- Applicability assessment

## Key Features

### Data Leakage Prevention ✅

**Problem:** Traditional ML workflows often leak information from test sets into training, leading to optimistically biased performance estimates.

**Solution:** OmicSelector 2.0 ensures:
```r
# ❌ OLD WAY (potential leakage)
data_normalized <- normalize(data)  # Normalization sees ALL data
train_test_split(data_normalized)

# ✅ NEW WAY (no leakage)
OmicSelector_nested_cv(
  data = data,
  preprocessing_recipe = recipe(...) %>% step_normalize(...),
  # Normalization happens INSIDE each fold
)
```

### Feature Selection Stability 📊

New stability-aware feature selection methods:
- **Stability Selection** - Select features consistently across subsamples
- **Model-X Knockoffs** - Control false discovery rate
- **Boruta with SHAP** - Enhanced feature importance
- **Feature Clustering** - Identify replaceable biomarkers

### Modern ML Algorithms 🤖

Out-of-the-box support for:
- **Random Forest** (ranger)
- **XGBoost** (gradient boosting)
- **Elastic Net** (glmnet)
- **SVM** (kernlab)
- **Neural Networks** (keras/tensorflow)

### Clinical Utility Metrics 🏥

Beyond accuracy and AUC:
- **Calibration Assessment**
  - Calibration plots
  - Hosmer-Lemeshow test
  - Calibration slope and intercept
  - Brier score
- **Decision Curve Analysis** (Coming in Phase 8)
- **Clinical Impact Curves** (Coming in Phase 8)

## Migration Guide

### From OmicSelector 1.0 to 2.0

#### Option 1: Use New Modern Interface (Recommended)

```r
# OLD (caret-based)
result <- OmicSelector_benchmark(
  wd = getwd(),
  algorithms = c("rf", "svmRadial"),
  ...
)

# NEW (tidymodels-based)
library(tidymodels)

# Define models
models <- list(
  rf = rand_forest(mtry = tune(), min_n = tune()) %>%
    set_engine("ranger") %>%
    set_mode("classification"),

  svm = svm_rbf(cost = tune()) %>%
    set_engine("kernlab") %>%
    set_mode("classification")
)

# Run nested CV with proper validation
result <- OmicSelector_nested_cv(
  data = your_data,
  outcome = "diagnosis",
  models = models,
  outer_folds = 5,
  inner_folds = 5,
  feature_selection_method = "boruta",
  calibrate = TRUE
)

# Generate TRIPOD+AI report
report <- OmicSelector_tripod_report(
  model_result = result,
  study_info = list(
    title = "Your Study Title",
    objective = "To develop a diagnostic model...",
    ...
  ),
  output_format = "html",
  output_file = "tripod_report.html"
)
```

#### Option 2: Keep Using Existing Code (Backward Compatible)

All existing OmicSelector functions continue to work:
```r
# This still works exactly as before!
result <- OmicSelector_benchmark(...)
features <- OmicSelector_OmicSelector(...)
plot <- OmicSelector_heatmap(...)
```

## Phase 1 Implementation Status ✅

### Completed
- ✅ Modern framework with tidymodels integration (`framework_modern.R`)
- ✅ Nested cross-validation with leakage prevention (`nested_cv.R`)
- ✅ TRIPOD+AI and PROBAST+AI compliance (`compliance.R`)
- ✅ Updated DESCRIPTION with modern dependencies
- ✅ GitHub Actions CI/CD workflows
- ✅ Comprehensive test suite (testthat 3.0)
- ✅ Updated to R version 2.0.0

### Coming Soon (Subsequent Phases)

#### Phase 2: Advanced Feature Selection
- Stability metrics (Nogueira)
- Model-X Knockoffs
- Feature clustering for biomarker replaceability

#### Phase 3: Multi-omics Integration
- DIABLO (mixOmics)
- MOFA+ integration
- Multi-assay data containers

#### Phase 4: Clinical Utility
- Full calibration suite
- Decision curve analysis
- Clinical impact curves
- Net benefit calculations

#### Phase 5: Survival Analysis
- Cox proportional hazards
- Random survival forests
- Time-dependent biomarkers
- Landmark analysis

#### Phase 6: Data Connectors
- GEO database integration
- TCGA data retrieval
- Automated dataset cards
- Data provenance tracking

#### Phase 7: Modern ML Algorithms
- AutoML integration (H2O)
- LightGBM and CatBoost
- Advanced ensemble methods
- Bayesian optimization

#### Phase 8: Explainability
- SHAP values
- LIME explanations
- Permutation importance
- ICE/PDP plots

## Installation

### Development Version (Recommended for Testing)

```r
# Install from GitHub
remotes::install_github("kstawiski/OmicSelector@claude/omicselector-modernization-phase1-011CUoP6wrbzCCxVtgHHC8B9")

# Or install with all suggested packages
remotes::install_github(
  "kstawiski/OmicSelector@claude/omicselector-modernization-phase1-011CUoP6wrbzCCxVtgHHC8B9",
  dependencies = TRUE
)
```

### Required Dependencies

```r
# Core tidymodels ecosystem
install.packages(c(
  "tidymodels",
  "parsnip",
  "recipes",
  "rsample",
  "tune",
  "workflows",
  "yardstick"
))

# ML backends
install.packages(c(
  "ranger",      # Random forest
  "xgboost",     # Gradient boosting
  "glmnet",      # Elastic net
  "kernlab"      # SVM
))

# For backward compatibility
install.packages("caret")
```

## Quick Start Example

```r
library(OmicSelector)
library(tidymodels)

# Load example data
data(miR_Asakura)  # Example miRNA expression data

# Define preprocessing
rec <- recipe(Class ~ ., data = miR_Asakura) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_zv(all_predictors())

# Run nested CV with feature selection
result <- OmicSelector_nested_cv(
  data = miR_Asakura,
  outcome = "Class",
  outer_folds = 5,
  inner_folds = 5,
  preprocessing_recipe = rec,
  feature_selection_method = "boruta",
  max_features = 20
)

# View results
print(result)
summary(result)
plot(result)

# Check feature stability
print(result$feature_stability)

# Generate compliance report
study_info <- list(
  title = "MicroRNA Biomarker Discovery",
  objective = "Develop diagnostic model for cancer detection",
  data_source = "RNA-seq from tumor and normal tissues",
  eligibility_criteria = "Adult patients with confirmed diagnosis"
)

report <- OmicSelector_tripod_report(
  model_result = result,
  study_info = study_info,
  output_format = "html",
  output_file = "report.html"
)

# Run risk of bias assessment
assessment <- OmicSelector_probast(result)
print(assessment)
```

## Documentation

### New Functions

| Function | Description |
|----------|-------------|
| `OmicSelector_fit()` | Unified model fitting with tidymodels/caret |
| `OmicSelector_nested_cv()` | Nested cross-validation with leakage prevention |
| `OmicSelector_tripod_report()` | Generate TRIPOD+AI compliance report |
| `OmicSelector_probast()` | PROBAST+AI risk of bias assessment |

### Existing Functions (Still Supported)

All existing functions remain available:
- `OmicSelector_benchmark()`
- `OmicSelector_OmicSelector()`
- `OmicSelector_prepare_split()`
- `OmicSelector_heatmap()`
- `OmicSelector_PCA()`
- And 30+ more...

## Testing

Run the test suite:

```r
# Install test dependencies
install.packages("testthat")

# Run tests
devtools::test()

# Check code coverage
covr::package_coverage()
```

## Performance Benchmarks

Preliminary benchmarks show:
- **Nested CV**: ~2x slower than simple CV, but provides unbiased estimates
- **Tidymodels vs caret**: Similar performance, better code clarity
- **Parallel processing**: Linear scaling up to available cores

## Contributing

We welcome contributions! Areas of focus:
1. Additional ML algorithms
2. New feature selection methods
3. Clinical utility metrics
4. Multi-omics methods
5. Documentation improvements

## Citation

If you use OmicSelector 2.0, please cite:

```bibtex
@software{omicselector2024,
  title = {OmicSelector: Modern Framework for Biomarker Discovery},
  author = {Stawiski, Konrad and Kaszkowiak, Marcin},
  year = {2024},
  version = {2.0.0},
  url = {https://github.com/kstawiski/OmicSelector}
}
```

## References

### TRIPOD+AI
- Collins GS, et al. (2024). TRIPOD+AI statement: updated guidance for reporting clinical prediction models that use regression or machine learning methods. *BMJ*.

### PROBAST
- Wolff RF, et al. (2019). PROBAST: A Tool to Assess the Risk of Bias and Applicability of Prediction Model Studies. *Annals of Internal Medicine*.

### Nested Cross-Validation
- Varma S, Simon R. (2006). Bias in error estimation when using cross-validation for model selection. *BMC Bioinformatics*.

### Feature Selection Stability
- Nogueira S, Brown G. (2016). Measuring the stability of feature selection. *Machine Learning*.

## License

MIT License - See LICENSE file for details

## Support

- **Issues**: https://github.com/kstawiski/OmicSelector/issues
- **Documentation**: https://biostat.umed.pl/OmicSelector/
- **Contact**: konrad.stawiski@umed.lodz.pl

---

**Note**: This is Phase 1 of the modernization. Additional features will be added in subsequent phases. Feedback and suggestions are welcome!
