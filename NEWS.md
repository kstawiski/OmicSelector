# OmicSelector 1.1.0 (Development)

## Major Changes

### Modern Machine Learning Framework
* **NEW**: Integration with tidymodels ecosystem for modern ML workflows
* **NEW**: `OmicSelector_fit()` - Unified interface supporting both tidymodels and caret
* **NEW**: `OmicSelector_nested_cv()` - Rigorous nested cross-validation with leakage prevention
* All preprocessing and feature selection now occurs within resampling folds by default
* Framework auto-detection with `method = "auto"`

### TRIPOD+AI and PROBAST+AI Compliance
* **NEW**: `OmicSelector_tripod_report()` - Automated TRIPOD+AI compliance reporting
* **NEW**: `OmicSelector_probast()` - PROBAST+AI risk of bias assessment
* **NEW**: `OmicSelector_model_card()` - Comprehensive model documentation
* Support for HTML, PDF, JSON, and DOCX output formats
* Automated checklist generation with completeness scoring

### Enhanced Validation and Reporting
* Nested cross-validation as default for unbiased performance estimation
* Automatic metadata tracking for reproducibility
* TRIPOD+AI compliant output includes all 27 required items
* PROBAST+AI assessment across 4 bias domains

## New Features

### Core Framework
* `OmicSelector_fit()` - Main model fitting function with framework routing
* `OmicSelector_nested_cv()` - Nested CV with proper validation
* Support for tidymodels parsnip model specifications
* Support for recipes preprocessing pipelines
* Automatic hyperparameter tuning within inner CV loops

### Compliance and Reporting
* TRIPOD+AI reporting with automated item detection
* PROBAST+AI risk of bias assessment
* Model cards following ML documentation best practices
* JSON/HTML export for compliance reports

### Infrastructure
* Comprehensive test suite with testthat
* Modern vignettes demonstrating workflows
* Updated DESCRIPTION with tidymodels dependencies
* S3 methods for print, summary, and plot

## Package Dependencies

### New Required Packages
* tidymodels (>= 1.1.0)
* rsample (>= 1.1.0)
* parsnip (>= 1.0.0)
* workflows (>= 1.1.0)
* tune (>= 1.1.0)
* yardstick (>= 1.1.0)
* dials (>= 1.1.0)
* jsonlite (>= 1.8.0)

### New Suggested Packages
* testthat (>= 3.0.0) - for testing
* mixOmics (>= 6.26.0) - for multi-omics integration
* survival (>= 3.5.0) - for survival analysis
* censored (>= 0.2.0) - for survival models
* DALEX (>= 2.4.0) - for explainability
* GEOquery (>= 2.70.0) - for GEO data access
* mlr3 (>= 0.16.0) - alternative ML framework
* h2o (>= 3.42.0) - for AutoML
* targets (>= 1.3.0) - for reproducible pipelines
* quarto - for modern documentation

## Documentation

### New Vignettes
* `modern_workflow.Rmd` - Complete guide to modernized OmicSelector
* Demonstrates nested CV, TRIPOD+AI, and best practices
* Examples of feature stability analysis
* Calibration and clinical utility assessment

### Updated Documentation
* All new functions have comprehensive roxygen2 documentation
* Examples provided for each major function
* References to TRIPOD+AI and PROBAST+AI guidelines

## Breaking Changes

**None** - This release maintains full backward compatibility with OmicSelector 1.0.0.
All existing functions continue to work as before. New functionality is additive.

## Improvements

### Code Quality
* Comprehensive test coverage for new functions
* Input validation for all public functions
* Consistent error messages and warnings
* S3 class system for result objects

### Performance
* Tidymodels backend optimized for large datasets
* Parallel processing support in nested CV
* Memory-efficient prediction storage options

### Reproducibility
* Comprehensive metadata tracking
* Seed management for all random operations
* Version tracking of R and package versions
* Complete preprocessing history

## Bug Fixes

* None in this release (new functionality)

## Known Issues

### Planned for Future Releases
* Full implementation of helper functions in `framework_modern.R`
* Multi-omics integration (DIABLO, MOFA) - coming in 1.2.0
* Advanced feature selection methods - coming in 1.2.0
* Survival analysis module - coming in 1.2.0
* Data connectors (GEO, TCGA) - coming in 1.2.0
* Model explainability (SHAP, LIME) - coming in 1.3.0
* AutoML integration - coming in 1.3.0

## Migration Guide

### From OmicSelector 1.0.0

No breaking changes! Your existing code will continue to work.

To take advantage of new features:

```r
# Old way (still works)
result <- OmicSelector_benchmark(...)

# New way (recommended for new projects)
library(tidymodels)

# Define preprocessing
rec <- recipe(outcome ~ ., data = data) %>%
  step_normalize(all_numeric_predictors())

# Define model
rf_spec <- rand_forest(mtry = tune()) %>%
  set_engine("ranger") %>%
  set_mode("classification")

# Run nested CV
result <- OmicSelector_nested_cv(
  data = data,
  outcome = "outcome",
  preprocessing_recipe = rec,
  models = list(rf = rf_spec),
  outer_folds = 5,
  inner_folds = 5
)

# Generate TRIPOD+AI report
tripod <- OmicSelector_tripod_report(
  model_result = result,
  study_info = list(title = "My Study", ...),
  output_file = "tripod_report.html"
)
```

## Acknowledgments

* TRIPOD+AI statement authors for transparent reporting guidelines
* PROBAST authors for bias assessment framework
* tidymodels team for excellent ML infrastructure
* OmicSelector users for feedback and suggestions

## References

1. Collins GS, Moons KGM, Dhiman P, et al. TRIPOD+AI statement: updated guidance for reporting clinical prediction models that use regression or machine learning methods. BMJ 2024;385:e078378. doi: 10.1136/bmj-2023-078378

2. Wolff RF, Moons KGM, Riley RD, et al. PROBAST: A Tool to Assess the Risk of Bias and Applicability of Prediction Model Studies. Ann Intern Med. 2019;170(1):51-58. doi:10.7326/M18-1376

3. Kuhn M, Wickham H. Tidymodels: a collection of packages for modeling and machine learning using tidyverse principles. https://www.tidymodels.org

---

# OmicSelector 1.0.0

Initial release of OmicSelector package for biomarker discovery from high-throughput omics experiments.

## Features

* 94 feature selection methods through 25 distinct approaches
* caret-based machine learning framework
* Keras/TensorFlow integration for deep learning
* Docker deployment with RStudio, Jupyter, VS Code
* Shiny web application
* ComBat for batch correction
* SMOTE for class balancing
* Comprehensive benchmarking functionality

See package documentation for details.
