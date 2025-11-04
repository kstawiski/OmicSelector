# OmicSelector Modernization Roadmap

This document tracks the progress of modernizing OmicSelector to 2025 standards.

## ✅ Phase 1: Core Architecture Setup (COMPLETED)

### Framework Layer
- [x] Create `R/framework_modern.R`
- [x] Implement `OmicSelector_fit()` with framework routing
- [x] Implement `OmicSelector_nested_cv()` with leakage prevention
- [x] Add tidymodels backend support
- [x] Maintain caret backward compatibility
- [x] Auto-detection of best framework
- [x] S3 methods (print, summary)

### Compliance Module
- [x] Create `R/compliance.R`
- [x] Implement `OmicSelector_tripod_report()` with 27 items
- [x] Implement `OmicSelector_probast()` for bias assessment
- [x] Implement `OmicSelector_model_card()` for documentation
- [x] Support HTML, JSON, PDF, DOCX export
- [x] Automatic metadata extraction
- [x] Completeness scoring
- [x] Recommendation generation

### Package Infrastructure
- [x] Update DESCRIPTION with tidymodels dependencies
- [x] Add suggested packages (mixOmics, survival, etc.)
- [x] Create test infrastructure (tests/testthat/)
- [x] Write comprehensive tests for nested CV
- [x] Write comprehensive tests for compliance
- [x] Create modern workflow vignette
- [x] Create NEWS.md
- [x] Create README_MODERNIZATION.md
- [x] Create GitHub Actions workflow for testing

## 🔄 Phase 2: Complete Core Implementation (IN PROGRESS)

### Helper Functions to Implement
The following helper functions in `R/framework_modern.R` need full implementation:

- [ ] `.fit_tidymodels()` - Complete tidymodels workflow
  - [ ] Integrate with workflows package
  - [ ] Proper tuning grid setup
  - [ ] Metric collection
  - [ ] Model fitting and evaluation

- [ ] `.fit_caret()` - Complete caret workflow
  - [ ] trainControl setup
  - [ ] Model training
  - [ ] Performance extraction

- [ ] `.perform_feature_selection()` - Feature selection methods
  - [ ] Boruta implementation
  - [ ] Stability selection
  - [ ] Variance-based filtering
  - [ ] Correlation-based filtering
  - [ ] Hook for custom methods

- [ ] `.tune_model_inner_cv()` - Inner CV tuning
  - [ ] Grid search
  - [ ] Random search
  - [ ] Bayesian optimization
  - [ ] Early stopping

- [ ] `.fit_final_model()` - Final model fitting
  - [ ] Parameter finalization
  - [ ] Model training on full outer fold
  - [ ] Importance extraction

- [ ] `.evaluate_predictions()` - Prediction evaluation
  - [ ] Classification metrics (AUC, accuracy, sensitivity, specificity)
  - [ ] Regression metrics (RMSE, MAE, R²)
  - [ ] Calibration metrics
  - [ ] Confidence intervals

- [ ] `.aggregate_nested_cv_results()` - Result aggregation
  - [ ] Mean and SD across folds
  - [ ] Confidence intervals
  - [ ] Best model selection
  - [ ] Stability metrics

- [ ] `.assess_calibration()` - Calibration assessment
  - [ ] Calibration plots
  - [ ] Brier score
  - [ ] Calibration slope/intercept
  - [ ] Hosmer-Lemeshow test

### Additional Core Functions
- [ ] `OmicSelector_calibrate()` - Standalone calibration
- [ ] `OmicSelector_decision_curve()` - Decision curve analysis
- [ ] `OmicSelector_clinical_impact()` - Clinical impact curves
- [ ] Plot methods for nested CV results
- [ ] Summary methods with detailed output

## 📊 Phase 3: Advanced Feature Selection

### Stability-Based Selection
- [ ] Create `R/feature_selection_modern.R`
- [ ] `OmicSelector_stable_features()` - Stability selection
- [ ] Nogueira stability metrics
- [ ] Bootstrap aggregating
- [ ] Subsampling approach
- [ ] Stability score calculation

### Statistical Methods
- [ ] Model-X knockoffs for FDR control
- [ ] Boruta with SHAP values
- [ ] Permutation importance
- [ ] Recursive feature elimination (updated)

### Biological Methods
- [ ] `OmicSelector_cluster_features()` - Feature clustering
- [ ] Pathway-based grouping
- [ ] Correlation-based clustering
- [ ] Representative selection
- [ ] Replaceability analysis

### Visualization
- [ ] Stability plots
- [ ] Feature importance plots
- [ ] Feature correlation networks
- [ ] Selection frequency heatmaps

## 🧬 Phase 4: Multi-omics Integration

### DIABLO Integration
- [ ] Create `R/multiomics.R`
- [ ] `OmicSelector_DIABLO()` - DIABLO wrapper
- [ ] Design matrix setup
- [ ] Multi-block PLS-DA
- [ ] Signature extraction
- [ ] Performance evaluation
- [ ] Integration plots (circos, arrow)

### MOFA+ Integration
- [ ] `OmicSelector_MOFA()` - MOFA+ wrapper
- [ ] Factor model fitting
- [ ] Variance decomposition
- [ ] Shared vs specific factors
- [ ] Factor interpretation

### Data Management
- [ ] `OmicSelector_MultiAssayExperiment()` - Data container
- [ ] Sample alignment checking
- [ ] Missing data handling
- [ ] Batch effect correction
- [ ] Quality control plots

### Validation
- [ ] Cross-omics validation
- [ ] Signature stability
- [ ] Platform transferability

## 🏥 Phase 5: Clinical Utility Metrics

### Calibration
- [ ] Create `R/clinical_utility.R`
- [ ] `OmicSelector_calibrate()` - Full implementation
- [ ] Platt scaling
- [ ] Isotonic regression
- [ ] Beta calibration
- [ ] Calibration plots (observed vs predicted)
- [ ] Calibration metrics:
  - [ ] Calibration slope
  - [ ] Calibration intercept
  - [ ] Brier score
  - [ ] ICI (Integrated Calibration Index)
  - [ ] E-statistics (Emax, E90, E50)

### Decision Analysis
- [ ] `OmicSelector_decision_curve()` - Full implementation
- [ ] Net benefit calculation
- [ ] Treat-all strategy
- [ ] Treat-none strategy
- [ ] Decision curve plotting
- [ ] Threshold optimization

### Clinical Impact
- [ ] `OmicSelector_clinical_impact()` - Full implementation
- [ ] Number needed to screen
- [ ] True/false positives by threshold
- [ ] Clinical impact curves
- [ ] Cost-effectiveness analysis (optional)

### Risk Stratification
- [ ] `OmicSelector_risk_groups()` - Risk grouping
- [ ] Cutpoint optimization
- [ ] Survival curves by risk group
- [ ] Reclassification tables
- [ ] Net reclassification improvement (NRI)
- [ ] Integrated discrimination improvement (IDI)

## ⏱️ Phase 6: Survival Analysis

### Core Survival Models
- [ ] Create `R/survival.R`
- [ ] `OmicSelector_survival()` - Survival model fitting
- [ ] Cox proportional hazards (censored package)
- [ ] Random survival forests
- [ ] XGBoost for survival (xgboost_cox)
- [ ] Deep learning survival (deepsurv)

### Survival Metrics
- [ ] C-index (Harrell's concordance)
- [ ] Time-dependent AUC
- [ ] Integrated Brier Score (IBS)
- [ ] Cumulative AUC

### Advanced Survival
- [ ] `OmicSelector_time_dependent()` - Time-varying effects
- [ ] Landmark analysis
- [ ] Competing risks
- [ ] Multi-state models

### Visualization
- [ ] Kaplan-Meier curves
- [ ] Forest plots
- [ ] Time-dependent ROC curves
- [ ] Calibration plots for survival

## 🔌 Phase 7: Data Connectors

### GEO Integration
- [ ] Create `R/data_connectors.R`
- [ ] `OmicSelector_GEO()` - GEO downloader
- [ ] GEOquery integration
- [ ] Automatic preprocessing
- [ ] Platform annotation
- [ ] Quality control

### TCGA Integration
- [ ] `OmicSelector_TCGA()` - TCGA accessor
- [ ] TCGAbiolinks integration
- [ ] Multi-omics download
- [ ] Clinical data integration
- [ ] Sample type filtering

### Dataset Documentation
- [ ] `OmicSelector_dataset_card()` - Dataset documentation
- [ ] Study description
- [ ] Sample characteristics
- [ ] Inclusion/exclusion criteria
- [ ] Preprocessing history
- [ ] Quality metrics
- [ ] Data provenance

### Format Converters
- [ ] `OmicSelector_to_MAE()` - Convert to MultiAssayExperiment
- [ ] `OmicSelector_to_SE()` - Convert to SummarizedExperiment
- [ ] `OmicSelector_from_MAE()` - Import from MAE
- [ ] `OmicSelector_from_SE()` - Import from SE

## 🤖 Phase 8: Modern ML Algorithms

### AutoML
- [ ] Create `R/ml_modern.R`
- [ ] `OmicSelector_automl()` - AutoML wrapper
- [ ] H2O AutoML integration
- [ ] Time budget management
- [ ] Algorithm selection
- [ ] Ensemble creation

### Ensemble Methods
- [ ] `OmicSelector_ensemble()` - Ensemble builder
- [ ] Stacking with meta-learner
- [ ] Voting (hard/soft)
- [ ] Bayesian model averaging
- [ ] Dynamic ensemble selection

### Deep Learning Updates
- [ ] Update Keras/TensorFlow integration
- [ ] Modern architectures
- [ ] Attention mechanisms
- [ ] Transformer models
- [ ] AutoML for DL (Auto-Keras)

### Specialized Algorithms
- [ ] LightGBM integration
- [ ] CatBoost integration
- [ ] NGBoost for uncertainty
- [ ] Quantile regression forests

## 🔍 Phase 9: Explainability

### Global Explanations
- [ ] Create `R/explainability.R`
- [ ] `OmicSelector_explain()` - Unified explainer
- [ ] SHAP values (fastshap/shapr)
- [ ] Permutation importance
- [ ] Partial dependence plots (PDP)
- [ ] Individual conditional expectation (ICE)

### Local Explanations
- [ ] LIME for individual predictions
- [ ] Local SHAP
- [ ] Counterfactual explanations
- [ ] Anchors

### Feature Importance
- [ ] `OmicSelector_importance()` - Aggregated importance
- [ ] Model-agnostic importance
- [ ] Grouped importance
- [ ] Conditional importance
- [ ] Confidence intervals

### Visualization
- [ ] SHAP summary plots
- [ ] SHAP force plots
- [ ] SHAP dependence plots
- [ ] Feature interaction plots
- [ ] Decision tree surrogate models

## 📦 Phase 10: Package Polish

### Documentation
- [ ] Complete all roxygen2 documentation
- [ ] Add examples to all functions
- [ ] Create additional vignettes:
  - [ ] Multi-omics integration
  - [ ] Survival analysis
  - [ ] Advanced feature selection
  - [ ] Clinical utility assessment
  - [ ] Explainability
  - [ ] Migration from 1.0
- [ ] Update pkgdown website
- [ ] Create cheat sheet
- [ ] Video tutorials

### Testing
- [ ] Achieve >80% code coverage
- [ ] Add tests for all helper functions
- [ ] Integration tests for workflows
- [ ] Performance benchmarks
- [ ] Edge case handling
- [ ] Memory leak testing

### Performance
- [ ] Optimize nested CV for large datasets
- [ ] Parallel processing enhancements
- [ ] Memory efficiency improvements
- [ ] Caching strategies
- [ ] Progress bars and logging

### User Experience
- [ ] Better error messages
- [ ] Input validation
- [ ] Informative warnings
- [ ] Progress tracking
- [ ] Verbose/quiet modes
- [ ] Dry-run option

## 🔧 Phase 11: Migration Support

### Migration Tools
- [ ] Create `R/migration.R`
- [ ] `OmicSelector_migrate()` - Code converter
- [ ] Parse old caret code
- [ ] Convert to tidymodels
- [ ] Generate migration report
- [ ] Side-by-side comparison

### Backward Compatibility
- [ ] Ensure all 1.0 functions work
- [ ] Add deprecation warnings (gentle)
- [ ] Provide migration paths
- [ ] Document breaking changes (if any)

### Migration Guide
- [ ] Detailed migration vignette
- [ ] Common patterns conversion
- [ ] Troubleshooting guide
- [ ] FAQ section

## 🎯 Priority Order

### Immediate (Weeks 1-2) - ✅ COMPLETED
1. ✅ Core framework setup
2. ✅ Nested CV implementation
3. ✅ TRIPOD+AI compliance
4. ✅ Basic testing
5. ✅ Initial documentation

### Short-term (Weeks 3-4) - 🔄 CURRENT
1. 🔄 Complete helper functions in framework_modern.R
2. ⏳ Implement calibration methods
3. ⏳ Add feature selection methods
4. ⏳ Comprehensive testing
5. ⏳ Polish documentation

### Medium-term (Weeks 5-8)
1. ⏳ Multi-omics integration
2. ⏳ Survival analysis module
3. ⏳ Clinical utility metrics
4. ⏳ Data connectors
5. ⏳ Advanced ML algorithms

### Long-term (Weeks 9-12)
1. ⏳ Explainability module
2. ⏳ AutoML integration
3. ⏳ Performance optimization
4. ⏳ Complete documentation
5. ⏳ Migration tools

## 📋 Quality Checklist

Before marking any phase as complete:

- [ ] All functions have roxygen2 documentation
- [ ] Examples provided and tested
- [ ] Unit tests written and passing
- [ ] Integration tests for workflows
- [ ] Code follows tidyverse style guide
- [ ] No data leakage in CV procedures
- [ ] Proper error handling
- [ ] Informative error messages
- [ ] Function works on edge cases
- [ ] Memory efficient
- [ ] Reproducible results (seeds)
- [ ] Vignette section completed
- [ ] NEWS.md updated

## 🐛 Known Issues to Address

1. Helper functions in framework_modern.R need implementation
2. Calibration functions not yet implemented
3. Decision curve analysis not yet implemented
4. Feature selection methods need expansion
5. Plot methods need implementation
6. Some edge cases not yet tested
7. Performance optimization needed for large datasets
8. Parallel processing not fully optimized

## 💡 Ideas for Future Enhancements

### Beyond Initial Roadmap
- Interactive web app for model exploration
- RESTful API for predictions
- Model registry and versioning
- Experiment tracking (MLflow integration)
- Automated report generation
- Containerized deployment
- Cloud platform integration
- Federated learning support
- Privacy-preserving ML
- Fairness metrics and bias detection

## 📚 Resources

### Key References
- TRIPOD+AI: BMJ 2024;385:e078378
- PROBAST: Ann Intern Med 2019;170(1):51-58
- Tidymodels: https://www.tidymodels.org
- Feature Stability: JMLR 2017;18(174):1-54
- DIABLO: Bioinformatics 2019;35(17):3055-3062
- MOFA+: Genome Biol 2020;21:111

### Learning Resources
- Tidy Modeling with R: https://www.tmwr.org
- Feature Engineering and Selection: http://www.feat.engineering
- Applied Predictive Modeling: Kuhn & Johnson
- Clinical Prediction Models: Steyerberg

---

Last updated: 2025-11-04

**Note**: This is a living document. Update as progress is made.
