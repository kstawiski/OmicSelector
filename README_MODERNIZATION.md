# OmicSelector Modernization (2025)

## 🚀 What's New

OmicSelector has been modernized to meet 2025 standards for biomarker discovery and machine learning best practices. The package now includes:

### ✨ Key Features

1. **Modern ML Framework**
   - Tidymodels integration with backward caret compatibility
   - Unified interface: `OmicSelector_fit()` and `OmicSelector_nested_cv()`
   - Support for 100+ algorithms through parsnip
   - Recipe-based preprocessing pipelines

2. **Rigorous Validation**
   - Nested cross-validation to prevent data leakage
   - All preprocessing inside resampling folds
   - Proper hyperparameter tuning workflow
   - Unbiased performance estimation

3. **TRIPOD+AI Compliance**
   - Automated TRIPOD+AI reporting
   - PROBAST+AI risk of bias assessment
   - Model cards for documentation
   - 27-item checklist with auto-completion

4. **Clinical Utility**
   - Calibration assessment (coming soon)
   - Decision curve analysis (coming soon)
   - Net benefit calculations (coming soon)
   - Clinical impact curves (coming soon)

## 📦 Installation

```r
# Install from GitHub
remotes::install_github("kstawiski/OmicSelector")

# Load package
library(OmicSelector)
library(tidymodels)
```

## 🎯 Quick Start

```r
library(OmicSelector)
library(tidymodels)

# 1. Prepare data
data <- your_omics_data  # rows = samples, columns = features + outcome

# 2. Define preprocessing
rec <- recipe(outcome ~ ., data = data) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_zv(all_predictors())

# 3. Define models
rf_spec <- rand_forest(mtry = tune(), trees = 1000) %>%
  set_engine("ranger") %>%
  set_mode("classification")

# 4. Run nested CV (prevents data leakage!)
result <- OmicSelector_nested_cv(
  data = data,
  outcome = "outcome",
  outer_folds = 5,
  inner_folds = 5,
  preprocessing_recipe = rec,
  models = list(random_forest = rf_spec),
  calibrate = TRUE,
  seed = 2024
)

# 5. View results
print(result)
summary(result)

# 6. Generate TRIPOD+AI report
tripod <- OmicSelector_tripod_report(
  model_result = result,
  study_info = list(
    title = "My Study",
    study_design = "cross-sectional",
    outcome_definition = "Disease status"
  ),
  output_file = "tripod_report.html"
)

# 7. Assess risk of bias
probast <- OmicSelector_probast(
  model_result = result,
  output_file = "probast_assessment.html"
)
```

## 🔬 Why This Matters

### Problem: Data Leakage

Traditional workflows often preprocess data before cross-validation:

```r
# ❌ BAD: Data leakage!
data_normalized <- scale(data)  # Uses info from ALL samples
cv_result <- cv_model(data_normalized)  # Optimistic bias!
```

This causes **optimistic bias** - your model looks better than it really is.

### Solution: Nested CV

OmicSelector ensures preprocessing happens inside folds:

```r
# ✅ GOOD: No leakage
result <- OmicSelector_nested_cv(
  data = data,  # Raw data
  preprocessing_recipe = recipe(...) %>% step_normalize(...),
  # Normalization applied separately in each fold
)
```

## 📊 What You Get

### Nested CV Results
- **Outer loop**: Unbiased performance estimates
- **Inner loop**: Hyperparameter tuning
- **Feature stability**: Which features are consistently selected
- **Model comparison**: Compare multiple algorithms fairly

### TRIPOD+AI Report
- 27-item checklist for transparent reporting
- Automated item detection from your analysis
- Recommendations for missing items
- HTML/PDF/JSON export

### PROBAST+AI Assessment
- Risk of bias across 4 domains:
  - Participants
  - Predictors
  - Outcome
  - Analysis
- Automatic detection of good practices (e.g., nested CV)
- Overall risk rating: Low/Moderate/High

## 🔄 Backward Compatibility

**All existing OmicSelector code continues to work!**

```r
# Old code still works
result <- OmicSelector_benchmark(...)

# But you can now also use modern features
result_modern <- OmicSelector_nested_cv(...)
```

You can gradually migrate to the new workflow.

## 📚 Documentation

### Vignettes
- `vignette("modern_workflow")` - Complete guide
- `vignette("multiomics_integration")` - Multi-omics (coming soon)
- `vignette("survival_analysis")` - Time-to-event (coming soon)

### Key Functions
- `OmicSelector_fit()` - Fit single model
- `OmicSelector_nested_cv()` - Rigorous nested CV
- `OmicSelector_tripod_report()` - TRIPOD+AI compliance
- `OmicSelector_probast()` - Bias assessment
- `OmicSelector_model_card()` - Model documentation

## 🗺️ Roadmap

### Phase 1: Core Architecture ✅ (Current)
- [x] Tidymodels integration
- [x] Nested cross-validation
- [x] TRIPOD+AI reporting
- [x] PROBAST+AI assessment
- [x] Comprehensive testing
- [x] Modern workflow vignette

### Phase 2: Feature Selection (Next)
- [ ] Stability selection
- [ ] Model-X knockoffs
- [ ] Feature clustering
- [ ] Stability metrics

### Phase 3: Multi-omics Integration
- [ ] DIABLO wrapper
- [ ] MOFA+ integration
- [ ] Multi-assay data containers
- [ ] Integrated signatures

### Phase 4: Clinical Utility
- [ ] Calibration methods
- [ ] Decision curve analysis
- [ ] Clinical impact curves
- [ ] Net benefit calculations

### Phase 5: Survival Analysis
- [ ] Cox models with tidymodels
- [ ] Random survival forests
- [ ] Time-dependent AUC
- [ ] Landmark analysis

### Phase 6: Data Connectors
- [ ] GEO database integration
- [ ] TCGA data access
- [ ] Dataset cards
- [ ] Automated preprocessing

### Phase 7: Advanced ML
- [ ] AutoML integration (H2O)
- [ ] Ensemble methods
- [ ] Deep learning updates
- [ ] GPU optimization

### Phase 8: Explainability
- [ ] SHAP values
- [ ] LIME explanations
- [ ] Permutation importance
- [ ] ICE/PDP plots

## 🧪 Testing

The package now includes comprehensive tests:

```r
# Run all tests
devtools::test()

# Run specific test file
testthat::test_file("tests/testthat/test-nested_cv.R")
```

Current test coverage: **framework and compliance modules**

## 🤝 Contributing

Contributions welcome! Please see CONTRIBUTING.md (coming soon).

Key areas for contribution:
- Helper function implementation
- Additional feature selection methods
- Multi-omics integration
- Clinical utility metrics
- Documentation improvements

## 📖 References

1. **TRIPOD+AI**: Collins GS, et al. BMJ 2024;385:e078378
2. **PROBAST**: Wolff RF, et al. Ann Intern Med 2019;170(1):51-58
3. **Tidymodels**: Kuhn M, Wickham H. https://www.tidymodels.org
4. **Feature Stability**: Nogueira S, et al. JMLR 2017;18(174):1-54

## 📄 License

MIT License - see LICENSE file

## 👥 Authors

- Konrad Stawiski (ORCID: 0000-0002-6550-3384)
- Marcin Kaszkowiak
- Contributors welcome!

## 📧 Contact

- Issues: https://github.com/kstawiski/OmicSelector/issues
- Email: konrad.stawiski@umed.lodz.pl
- Website: https://biostat.umed.pl/OmicSelector/

---

**Note**: This is a development version. Some features are still being implemented.
See NEWS.md for detailed changelog and roadmap.
