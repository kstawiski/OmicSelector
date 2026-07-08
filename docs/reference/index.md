# Package index

## Core Pipeline

Main classes for building zero-leakage ML pipelines

- [`OmicPipeline`](https://kstawiski.github.io/OmicSelector/reference/OmicPipeline.md)
  : OmicPipeline: Zero-Leakage Feature Selection Pipeline
- [`omic_pipeline()`](https://kstawiski.github.io/OmicSelector/reference/omic_pipeline.md)
  : Quick pipeline creation from data
- [`BenchmarkService`](https://kstawiski.github.io/OmicSelector/reference/BenchmarkService.md)
  : BenchmarkService: Nested Cross-Validation with Zero Leakage
- [`omic_benchmark()`](https://kstawiski.github.io/OmicSelector/reference/omic_benchmark.md)
  : Create benchmark service from OmicPipeline

## Signature Selection

Multi-objective biomarker signature selection

- [`select_best_signature()`](https://kstawiski.github.io/OmicSelector/reference/select_best_signature.md)
  : Select Best Biomarker Signature from Nested CV Results
- [`get_consensus_features()`](https://kstawiski.github.io/OmicSelector/reference/get_consensus_features.md)
  : Get Consensus Features from Best Signature
- [`get_selected_features_per_fold()`](https://kstawiski.github.io/OmicSelector/reference/get_selected_features_per_fold.md)
  : Get Selected Features Per Fold

## Stability Analysis

Feature selection stability metrics

- [`compute_nogueira_stability()`](https://kstawiski.github.io/OmicSelector/reference/compute_nogueira_stability.md)
  : Compute Nogueira Stability Index
- [`compute_stability_from_resample()`](https://kstawiski.github.io/OmicSelector/reference/compute_stability_from_resample.md)
  : Compute Stability from ResampleResult
- [`extract_features_from_resample()`](https://kstawiski.github.io/OmicSelector/reference/extract_features_from_resample.md)
  : Extract Features from All Folds in ResampleResult
- [`extract_selected_features()`](https://kstawiski.github.io/OmicSelector/reference/extract_selected_features.md)
  : Extract Selected Features from Trained GraphLearner

## Batch Correction

FrozenComBat for proper batch effect handling

- [`FrozenComBat`](https://kstawiski.github.io/OmicSelector/reference/FrozenComBat.md)
  : FrozenComBat R6 Class
- [`frozen_combat_correct()`](https://kstawiski.github.io/OmicSelector/reference/frozen_combat_correct.md)
  : Convenience function for frozen ComBat correction
- [`create_frozen_combat_pipeop()`](https://kstawiski.github.io/OmicSelector/reference/create_frozen_combat_pipeop.md)
  : Create a Frozen ComBat PipeOp for mlr3pipelines

## Cross-Platform Transfer

Cross-platform model deployment and domain adaptation

- [`cross_platform_transfer()`](https://kstawiski.github.io/OmicSelector/reference/cross_platform_transfer.md)
  : Convenience Wrapper for Cross-Platform Transfer

## Within-Sample CoDA Methods

Paper 1 v2.2 image encodings, neural learners, and perturbation
benchmarking

- [`encode_simple_grid()`](https://kstawiski.github.io/OmicSelector/reference/encode_simple_grid.md)
  : Simple Row-Major Grid Encoding
- [`encode_corr_grid()`](https://kstawiski.github.io/OmicSelector/reference/encode_corr_grid.md)
  : Encode with Correlation-Ordered Grid
- [`encode_deepinsight()`](https://kstawiski.github.io/OmicSelector/reference/encode_deepinsight.md)
  : Encode with DeepInsight Layout
- [`encode_ratio_image()`](https://kstawiski.github.io/OmicSelector/reference/encode_ratio_image.md)
  : Create Pairwise Log-Ratio Image
- [`train_ratio_cnn()`](https://kstawiski.github.io/OmicSelector/reference/train_ratio_cnn.md)
  : Train Ratio Image CNN
- [`train_ratio_cnn_multiseed()`](https://kstawiski.github.io/OmicSelector/reference/train_ratio_cnn_multiseed.md)
  : Train Ratio CNN Across Three Default Seeds
- [`clr_transform()`](https://kstawiski.github.io/OmicSelector/reference/clr_transform.md)
  : Centered Log-Ratio Transformation
- [`train_clr_mlp()`](https://kstawiski.github.io/OmicSelector/reference/train_clr_mlp.md)
  : Train a CLR + MLP Classifier
- [`predict_clr_mlp()`](https://kstawiski.github.io/OmicSelector/reference/predict_clr_mlp.md)
  : Predict with a CLR + MLP Classifier
- [`train_codacore()`](https://kstawiski.github.io/OmicSelector/reference/train_codacore.md)
  : Train a CoDaCoRe Classifier
- [`predict_codacore()`](https://kstawiski.github.io/OmicSelector/reference/predict_codacore.md)
  : Predict with a CoDaCoRe Classifier
- [`ws_perturbation_benchmark()`](https://kstawiski.github.io/OmicSelector/reference/ws_perturbation_benchmark.md)
  : Within-Sample Perturbation Benchmark

## Calibration

Probability calibration and diagnostics

- [`fit_platt_scaling()`](https://kstawiski.github.io/OmicSelector/reference/fit_platt_scaling.md)
  : Platt Scaling (Logistic Calibration)
- [`fit_isotonic_calibration()`](https://kstawiski.github.io/OmicSelector/reference/fit_isotonic_calibration.md)
  : Isotonic Regression Calibration
- [`fit_temperature_scaling()`](https://kstawiski.github.io/OmicSelector/reference/fit_temperature_scaling.md)
  : Temperature Scaling
- [`compute_ece()`](https://kstawiski.github.io/OmicSelector/reference/compute_ece.md)
  : Compute Expected Calibration Error (ECE)
- [`decompose_brier()`](https://kstawiski.github.io/OmicSelector/reference/decompose_brier.md)
  : Decompose Brier Score
- [`calibration_summary()`](https://kstawiski.github.io/OmicSelector/reference/calibration_summary.md)
  : Calibration Summary for Model Results
- [`reliability_diagram_data()`](https://kstawiski.github.io/OmicSelector/reference/reliability_diagram_data.md)
  : Create Reliability Diagram Data

## Interpretability

Model interpretation and feature importance

- [`create_explainer()`](https://kstawiski.github.io/OmicSelector/reference/create_explainer.md)
  : Create Model Explainer
- [`shap_values()`](https://kstawiski.github.io/OmicSelector/reference/shap_values.md)
  : Compute SHAP-like Values
- [`feature_importance()`](https://kstawiski.github.io/OmicSelector/reference/feature_importance.md)
  : Compute Permutation Feature Importance
- [`partial_dependence()`](https://kstawiski.github.io/OmicSelector/reference/partial_dependence.md)
  : Compute Partial Dependence
- [`check_feature_correlations()`](https://kstawiski.github.io/OmicSelector/reference/check_feature_correlations.md)
  : Check Feature Correlations

## Phase 5: GOF Filters

Goodness-of-fit filters for sparse/zero-inflated omics data

- [`FilterGOF_KS`](https://kstawiski.github.io/OmicSelector/reference/FilterGOF_KS.md)
  : Kolmogorov-Smirnov GOF Filter
- [`FilterHurdle`](https://kstawiski.github.io/OmicSelector/reference/FilterHurdle.md)
  : Hurdle Filter for Zero-Inflated Data
- [`FilterZeroProp`](https://kstawiski.github.io/OmicSelector/reference/FilterZeroProp.md)
  : Zero-Proportion Filter
- [`make_gof_filter()`](https://kstawiski.github.io/OmicSelector/reference/make_gof_filter.md)
  : Create GOF Filter
- [`compare_gof_filters()`](https://kstawiski.github.io/OmicSelector/reference/compare_gof_filters.md)
  : Compare GOF Filters on Task
- [`register_gof_filters()`](https://kstawiski.github.io/OmicSelector/reference/register_gof_filters.md)
  : Register GOF Filters in mlr3

## Phase 5: Bayesian Tuning

Bayesian hyperparameter optimization with mlr3mbo

- [`make_autotuner_glmnet()`](https://kstawiski.github.io/OmicSelector/reference/make_autotuner_glmnet.md)
  : Create Bayesian-Optimized AutoTuner for glmnet
- [`make_autotuner_xgboost()`](https://kstawiski.github.io/OmicSelector/reference/make_autotuner_xgboost.md)
  : Create Bayesian-Optimized AutoTuner for XGBoost
- [`make_autotuner_ranger()`](https://kstawiski.github.io/OmicSelector/reference/make_autotuner_ranger.md)
  : Create Bayesian-Optimized AutoTuner for Random Forest
- [`make_autotuner_lightgbm()`](https://kstawiski.github.io/OmicSelector/reference/make_autotuner_lightgbm.md)
  : Create Bayesian-Optimized AutoTuner for LightGBM
- [`get_optimal_params()`](https://kstawiski.github.io/OmicSelector/reference/get_optimal_params.md)
  : Get Optimal Hyperparameters from AutoTuner
- [`run_bayesian_benchmark()`](https://kstawiski.github.io/OmicSelector/reference/run_bayesian_benchmark.md)
  : Run Bayesian Optimization Benchmark

## Phase 5: AutoXAI

DALEX-based interpretability with correlation warnings

- [`xai_pipeline()`](https://kstawiski.github.io/OmicSelector/reference/xai_pipeline.md)
  : Run Complete XAI Pipeline
- [`xai_explainer_mlr3()`](https://kstawiski.github.io/OmicSelector/reference/xai_explainer_mlr3.md)
  : Create DALEX Explainer from mlr3 Learner
- [`xai_importance()`](https://kstawiski.github.io/OmicSelector/reference/xai_importance.md)
  : Compute Permutation Feature Importance
- [`xai_correlations()`](https://kstawiski.github.io/OmicSelector/reference/xai_correlations.md)
  : Compute Correlation Diagnostics for Features
- [`xai_pdp()`](https://kstawiski.github.io/OmicSelector/reference/xai_pdp.md)
  : Compute Partial Dependence Plots
- [`xai_shap()`](https://kstawiski.github.io/OmicSelector/reference/xai_shap.md)
  : Compute SHAP Values for Observations
- [`plot_xai_importance()`](https://kstawiski.github.io/OmicSelector/reference/plot_xai_importance.md)
  : Plot XAI Feature Importance
- [`print_xai_summary()`](https://kstawiski.github.io/OmicSelector/reference/print_xai_summary.md)
  : Print XAI Summary

## Phase 5: Stability Ensemble

Bootstrap-based feature selection stability

- [`create_stability_ensemble()`](https://kstawiski.github.io/OmicSelector/reference/create_stability_ensemble.md)
  : Create a Stability Ensemble with Presets
- [`StabilityEnsemble`](https://kstawiski.github.io/OmicSelector/reference/StabilityEnsemble.md)
  : Stability-Based Ensemble Selection
- [`SequentialSelector`](https://kstawiski.github.io/OmicSelector/reference/SequentialSelector.md)
  : Hybrid Sequential Feature Selection (HSFS)

## Phase 5: Synthetic Data

Data augmentation and synthetic data generation

- [`smote_augment()`](https://kstawiski.github.io/OmicSelector/reference/smote_augment.md)
  : SMOTE Augmentation for Omics Data
- [`noise_augment()`](https://kstawiski.github.io/OmicSelector/reference/noise_augment.md)
  : Gaussian Noise Augmentation
- [`validate_synthetic()`](https://kstawiski.github.io/OmicSelector/reference/validate_synthetic.md)
  : Validate Synthetic Data Quality
- [`balance_classes()`](https://kstawiski.github.io/OmicSelector/reference/balance_classes.md)
  : Create Balanced Training Set
- [`tabddpm_generate()`](https://kstawiski.github.io/OmicSelector/reference/tabddpm_generate.md)
  : TabDDPM Synthetic Data Generator

## Phase 5: Deep Learning

Deep learning models for omics (requires torch)

- [`create_mlp_learner()`](https://kstawiski.github.io/OmicSelector/reference/create_mlp_learner.md)
  : Create MLP Learner via mlr3torch
- [`make_mlp_learner()`](https://kstawiski.github.io/OmicSelector/reference/make_mlp_learner.md)
  : Make MLP Learner (Alias)
- [`make_fttransformer_learner()`](https://kstawiski.github.io/OmicSelector/reference/make_fttransformer_learner.md)
  [`make_tabtransformer_learner()`](https://kstawiski.github.io/OmicSelector/reference/make_fttransformer_learner.md)
  : Create an FT-Transformer Learner
- [`make_gnn_learner()`](https://kstawiski.github.io/OmicSelector/reference/make_gnn_learner.md)
  : Create GNN Learner for Pathway-Aware Classification
- [`build_correlation_adjacency()`](https://kstawiski.github.io/OmicSelector/reference/build_correlation_adjacency.md)
  : Create Correlation-Based Adjacency for GNN
- [`run_dl_benchmark()`](https://kstawiski.github.io/OmicSelector/reference/run_dl_benchmark.md)
  : Deep Learning Benchmark
- [`check_dl_availability()`](https://kstawiski.github.io/OmicSelector/reference/check_dl_availability.md)
  : Check Deep Learning Availability

## Multi-Omics

Multi-omics late integration

- [`merge_omics_data()`](https://kstawiski.github.io/OmicSelector/reference/merge_omics_data.md)
  : Merge Multi-Omics Data for Analysis
- [`get_modality_info()`](https://kstawiski.github.io/OmicSelector/reference/get_modality_info.md)
  : Get Modality Information
- [`validate_omics_input()`](https://kstawiski.github.io/OmicSelector/reference/validate_omics_input.md)
  : Validate Multi-Omics Input
- [`stack_omics()`](https://kstawiski.github.io/OmicSelector/reference/stack_omics.md)
  : Create Multi-Omics Stacked Ensemble (Convenience Function)

## Model Export & Reporting

Model deployment and TRIPOD reporting

- [`generate_tripod_report()`](https://kstawiski.github.io/OmicSelector/reference/generate_tripod_report.md)
  : Generate TRIPOD+AI Report
- [`create_report_data()`](https://kstawiski.github.io/OmicSelector/reference/create_report_data.md)
  : Create Report Data Schema

## Utilities

Helper functions for parallel processing and caching

- [`setup_parallel()`](https://kstawiski.github.io/OmicSelector/reference/setup_parallel.md)
  : Configure Parallelization for OmicSelector
- [`with_parallel()`](https://kstawiski.github.io/OmicSelector/reference/with_parallel.md)
  : With Parallel Scope
- [`get_parallel_status()`](https://kstawiski.github.io/OmicSelector/reference/get_parallel_status.md)
  : Get Current Parallelization Status
- [`reset_parallel()`](https://kstawiski.github.io/OmicSelector/reference/reset_parallel.md)
  : Reset Parallelization to Sequential
- [`create_omic_cache()`](https://kstawiski.github.io/OmicSelector/reference/create_omic_cache.md)
  : Create Split-Aware Cache
- [`cached_filter()`](https://kstawiski.github.io/OmicSelector/reference/cached_filter.md)
  : Cached Filter Computation
- [`cache_stats()`](https://kstawiski.github.io/OmicSelector/reference/cache_stats.md)
  : Get Cache Statistics
- [`clear_cache()`](https://kstawiski.github.io/OmicSelector/reference/clear_cache.md)
  : Clear OmicSelector Cache

## Legacy Functions

Functions from OmicSelector 1.x (use with caution)

- [`OmicSelector_correlation_plot()`](https://kstawiski.github.io/OmicSelector/reference/OmicSelector_correlation_plot.md)
  : OmicSelector_correlation_plot
- [`OmicSelector_profileplot()`](https://kstawiski.github.io/OmicSelector/reference/OmicSelector_profileplot.md)
  : OmicSelector_profileplot
- [`OmicSelector_propensity_score_matching()`](https://kstawiski.github.io/OmicSelector/reference/OmicSelector_propensity_score_matching.md)
  : OmicSelector_propensity_score_matching
- [`OmicSelector_vulcano_plot()`](https://kstawiski.github.io/OmicSelector/reference/OmicSelector_vulcano_plot.md)
  : OmicSelector_vulcano_plot

## Internal Functions

Internal functions and helpers (exported for advanced use)

- [`AutoXAI`](https://kstawiski.github.io/OmicSelector/reference/AutoXAI.md)
  : Auto XAI: Automatic Explainability for Biomarker Models
- [`BayesianTuner`](https://kstawiski.github.io/OmicSelector/reference/BayesianTuner.md)
  : Bayesian Hyperparameter Optimization for Omics
- [`BenchmarkService`](https://kstawiski.github.io/OmicSelector/reference/BenchmarkService.md)
  : BenchmarkService: Nested Cross-Validation with Zero Leakage
- [`CrossPlatformAdapter`](https://kstawiski.github.io/OmicSelector/reference/CrossPlatformAdapter.md)
  : CrossPlatformAdapter R6 Class
- [`DeepLearners`](https://kstawiski.github.io/OmicSelector/reference/DeepLearners.md)
  : Deep Learning Learners for Omics Data
- [`FilterCoDA_CoDaCoRe`](https://kstawiski.github.io/OmicSelector/reference/FilterCoDA_CoDaCoRe.md)
  : CoDaCoRe Feature Filter
- [`FilterCoDA_LogcontrastLasso`](https://kstawiski.github.io/OmicSelector/reference/FilterCoDA_LogcontrastLasso.md)
  : Log-Contrast Lasso Filter
- [`FilterCoDA_PLRVariance`](https://kstawiski.github.io/OmicSelector/reference/FilterCoDA_PLRVariance.md)
  : Pairwise Log-Ratio Variance Filter
- [`FilterCoDA_Selbal`](https://kstawiski.github.io/OmicSelector/reference/FilterCoDA_Selbal.md)
  : Selbal-Style Forward Balance Filter
- [`FilterCoDA_StabilityLogratio`](https://kstawiski.github.io/OmicSelector/reference/FilterCoDA_StabilityLogratio.md)
  : Stability Log-Ratio Filter
- [`FilterGOF`](https://kstawiski.github.io/OmicSelector/reference/FilterGOF.md)
  : Goodness-of-Fit Filters for Sparse Omics Data
- [`FilterGOF_KS`](https://kstawiski.github.io/OmicSelector/reference/FilterGOF_KS.md)
  : Kolmogorov-Smirnov GOF Filter
- [`FilterHurdle`](https://kstawiski.github.io/OmicSelector/reference/FilterHurdle.md)
  : Hurdle Filter for Zero-Inflated Data
- [`FilterZeroProp`](https://kstawiski.github.io/OmicSelector/reference/FilterZeroProp.md)
  : Zero-Proportion Filter
- [`FrozenComBat`](https://kstawiski.github.io/OmicSelector/reference/FrozenComBat.md)
  : FrozenComBat R6 Class
- [`GoldStandard`](https://kstawiski.github.io/OmicSelector/reference/GoldStandard.md)
  : Gold Standard Synthetic Dataset for Leakage Detection
- [`OmicModalitySpec`](https://kstawiski.github.io/OmicSelector/reference/OmicModalitySpec.md)
  : OmicModalitySpec R6 Class
- [`OmicPipeline`](https://kstawiski.github.io/OmicSelector/reference/OmicPipeline.md)
  : OmicPipeline: Zero-Leakage Feature Selection Pipeline
- [`OmicSelector-data`](https://kstawiski.github.io/OmicSelector/reference/OmicSelector-data.md)
  [`detectable_in_serum`](https://kstawiski.github.io/OmicSelector/reference/OmicSelector-data.md)
  [`miRNAselector_tutorial_balanced_benchmark`](https://kstawiski.github.io/OmicSelector/reference/OmicSelector-data.md)
  [`miRNAselector_tutorial_balanced_dataset`](https://kstawiski.github.io/OmicSelector/reference/OmicSelector-data.md)
  [`miRNAselector_tutorial_balanced_mixed`](https://kstawiski.github.io/OmicSelector/reference/OmicSelector-data.md)
  [`original_TCGA_data`](https://kstawiski.github.io/OmicSelector/reference/OmicSelector-data.md)
  [`orginal_TCGA_data`](https://kstawiski.github.io/OmicSelector/reference/OmicSelector-data.md)
  : Datasets included with OmicSelector
- [`OmicSelector-package`](https://kstawiski.github.io/OmicSelector/reference/OmicSelector-package.md)
  [`OmicSelector`](https://kstawiski.github.io/OmicSelector/reference/OmicSelector-package.md)
  : OmicSelector: Zero-Leakage Biomarker Discovery Toolkit
- [`OmicSelector_correlation_plot()`](https://kstawiski.github.io/OmicSelector/reference/OmicSelector_correlation_plot.md)
  : OmicSelector_correlation_plot
- [`OmicSelector_profileplot()`](https://kstawiski.github.io/OmicSelector/reference/OmicSelector_profileplot.md)
  : OmicSelector_profileplot
- [`OmicSelector_propensity_score_matching()`](https://kstawiski.github.io/OmicSelector/reference/OmicSelector_propensity_score_matching.md)
  : OmicSelector_propensity_score_matching
- [`OmicSelector_vulcano_plot()`](https://kstawiski.github.io/OmicSelector/reference/OmicSelector_vulcano_plot.md)
  : OmicSelector_vulcano_plot
- [`OmicStackedEnsemble`](https://kstawiski.github.io/OmicSelector/reference/OmicStackedEnsemble.md)
  : OmicStackedEnsemble R6 Class
- [`OmicWeightedEnsemble`](https://kstawiski.github.io/OmicSelector/reference/OmicWeightedEnsemble.md)
  : Simple Weighted Averaging Ensemble
- [`SequentialSelector`](https://kstawiski.github.io/OmicSelector/reference/SequentialSelector.md)
  : Hybrid Sequential Feature Selection (HSFS)
- [`StabilityEnsemble`](https://kstawiski.github.io/OmicSelector/reference/StabilityEnsemble.md)
  : Stability-Based Ensemble Selection
- [`SyntheticData`](https://kstawiski.github.io/OmicSelector/reference/SyntheticData.md)
  : Synthetic Data Generation for Omics
- [`TabularDL`](https://kstawiski.github.io/OmicSelector/reference/TabularDL.md)
  : Modern Tabular Learners for OmicSelector
- [`apply_compositional_mahalanobis()`](https://kstawiski.github.io/OmicSelector/reference/apply_compositional_mahalanobis.md)
  : Apply compositional Mahalanobis distance to new samples
- [`apply_frozen_combat_cv()`](https://kstawiski.github.io/OmicSelector/reference/apply_frozen_combat_cv.md)
  : Apply Frozen ComBat Within Cross-Validation Folds
- [`apply_frozen_quantile()`](https://kstawiski.github.io/OmicSelector/reference/apply_frozen_quantile.md)
  : Apply frozen quantile calibration to new samples
- [`apply_frozen_ruv()`](https://kstawiski.github.io/OmicSelector/reference/apply_frozen_ruv.md)
  : Apply frozen RUV correction to new samples
- [`apply_hemolysis_prefilter()`](https://kstawiski.github.io/OmicSelector/reference/apply_hemolysis_prefilter.md)
  : Apply a Hemolysis Pre-Filter
- [`apply_hemolysis_rr()`](https://kstawiski.github.io/OmicSelector/reference/apply_hemolysis_rr.md)
  : Apply frozen hemolysis-RR correction to new samples
- [`apply_isolation_forest_logratio()`](https://kstawiski.github.io/OmicSelector/reference/apply_isolation_forest_logratio.md)
  : Apply isolation forest to new samples and flag anomalies
- [`apply_logistic_normal_eb()`](https://kstawiski.github.io/OmicSelector/reference/apply_logistic_normal_eb.md)
  : Apply frozen logistic-normal EB shrinkage to new samples
- [`apply_mirna_aliases()`](https://kstawiski.github.io/OmicSelector/reference/apply_mirna_aliases.md)
  : Rename miRNA features in a matrix or named vector
- [`apply_robust_pca_residual()`](https://kstawiski.github.io/OmicSelector/reference/apply_robust_pca_residual.md)
  : Apply MAD-scaled SVD residual filter to new samples
- [`autoencoder`](https://kstawiski.github.io/OmicSelector/reference/autoencoder.md)
  : Autoencoder Utilities (torch)
- [`autoencoder_encode()`](https://kstawiski.github.io/OmicSelector/reference/autoencoder_encode.md)
  : Encode Features with Autoencoder
- [`autoencoder_fit()`](https://kstawiski.github.io/OmicSelector/reference/autoencoder_fit.md)
  : Fit Autoencoder
- [`autoencoder_load()`](https://kstawiski.github.io/OmicSelector/reference/autoencoder_load.md)
  : Load Autoencoder State
- [`autoencoder_save()`](https://kstawiski.github.io/OmicSelector/reference/autoencoder_save.md)
  : Save Autoencoder State
- [`balance_classes()`](https://kstawiski.github.io/OmicSelector/reference/balance_classes.md)
  : Create Balanced Training Set
- [`bias-audit`](https://kstawiski.github.io/OmicSelector/reference/bias-audit.md)
  : Cohort-Provenance Bias Audit
- [`build_correlation_adjacency()`](https://kstawiski.github.io/OmicSelector/reference/build_correlation_adjacency.md)
  : Create Correlation-Based Adjacency for GNN
- [`cache`](https://kstawiski.github.io/OmicSelector/reference/cache.md)
  : Split-Aware Caching for OmicSelector
- [`cache_stats()`](https://kstawiski.github.io/OmicSelector/reference/cache_stats.md)
  : Get Cache Statistics
- [`cached_filter()`](https://kstawiski.github.io/OmicSelector/reference/cached_filter.md)
  : Cached Filter Computation
- [`calibration`](https://kstawiski.github.io/OmicSelector/reference/calibration.md)
  : Calibration Metrics for OmicSelector
- [`calibration_summary()`](https://kstawiski.github.io/OmicSelector/reference/calibration_summary.md)
  : Calibration Summary for Model Results
- [`check_batch_correction_leakage()`](https://kstawiski.github.io/OmicSelector/reference/check_batch_correction_leakage.md)
  : Check for Batch Correction Leakage
- [`check_dl_availability()`](https://kstawiski.github.io/OmicSelector/reference/check_dl_availability.md)
  : Check Deep Learning Availability
- [`check_feature_correlations()`](https://kstawiski.github.io/OmicSelector/reference/check_feature_correlations.md)
  : Check Feature Correlations
- [`check_no_premature_imputation()`](https://kstawiski.github.io/OmicSelector/reference/check_no_premature_imputation.md)
  : Detect Global Pre-Loop Imputation (Leakage Check)
- [`check_null_benchmark_draws()`](https://kstawiski.github.io/OmicSelector/reference/check_null_benchmark_draws.md)
  : Validate Random-Panel Null Benchmark
- [`clear_cache()`](https://kstawiski.github.io/OmicSelector/reference/clear_cache.md)
  : Clear OmicSelector Cache
- [`clr-mlp`](https://kstawiski.github.io/OmicSelector/reference/clr-mlp.md)
  : CLR + MLP Models for Within-Sample Classification
- [`clr_transform()`](https://kstawiski.github.io/OmicSelector/reference/clr_transform.md)
  : Centered Log-Ratio Transformation
- [`coda-feature-selection`](https://kstawiski.github.io/OmicSelector/reference/coda-feature-selection.md)
  : CoDA-Aware Feature Selection for Biomarker Panels
- [`codaFS_codacore_wrapper()`](https://kstawiski.github.io/OmicSelector/reference/codaFS_codacore_wrapper.md)
  : CoDaCoRe-Derived Feature Selection
- [`codaFS_logcontrast_lasso()`](https://kstawiski.github.io/OmicSelector/reference/codaFS_logcontrast_lasso.md)
  : Log-Contrast Lasso on CLR Features
- [`codaFS_plr_variance()`](https://kstawiski.github.io/OmicSelector/reference/codaFS_plr_variance.md)
  : Pairwise Log-Ratio Variance Feature Filter
- [`codaFS_selbal_wrapper()`](https://kstawiski.github.io/OmicSelector/reference/codaFS_selbal_wrapper.md)
  : Selbal-Style Forward Balance Feature Selection
- [`codaFS_stability_logratio()`](https://kstawiski.github.io/OmicSelector/reference/codaFS_stability_logratio.md)
  : Stability Selection on Pairwise Log-Ratios
- [`codacore-interface`](https://kstawiski.github.io/OmicSelector/reference/codacore-interface.md)
  : CoDaCoRe Interfaces for Sparse Log-Contrast Classification
- [`compare_gof_filters()`](https://kstawiski.github.io/OmicSelector/reference/compare_gof_filters.md)
  : Compare GOF Filters on Task
- [`compute_domain_shift()`](https://kstawiski.github.io/OmicSelector/reference/compute_domain_shift.md)
  : Quantify Domain Shift Between Source and Target Platforms
- [`compute_ece()`](https://kstawiski.github.io/OmicSelector/reference/compute_ece.md)
  : Compute Expected Calibration Error (ECE)
- [`compute_nogueira_stability()`](https://kstawiski.github.io/OmicSelector/reference/compute_nogueira_stability.md)
  : Compute Nogueira Stability Index
- [`compute_shap_with_warnings()`](https://kstawiski.github.io/OmicSelector/reference/compute_shap_with_warnings.md)
  : Compute SHAP Values with Correlation Warnings
- [`compute_stability_from_resample()`](https://kstawiski.github.io/OmicSelector/reference/compute_stability_from_resample.md)
  : Compute Stability from ResampleResult
- [`create_autoencoder_pipeop()`](https://kstawiski.github.io/OmicSelector/reference/create_autoencoder_pipeop.md)
  : Create Autoencoder PipeOp
- [`create_explainer()`](https://kstawiski.github.io/OmicSelector/reference/create_explainer.md)
  : Create Model Explainer
- [`create_frozen_combat_pipeop()`](https://kstawiski.github.io/OmicSelector/reference/create_frozen_combat_pipeop.md)
  : Create a Frozen ComBat PipeOp for mlr3pipelines
- [`create_hsfs_selector()`](https://kstawiski.github.io/OmicSelector/reference/create_hsfs_selector.md)
  : Create a Hybrid Sequential Feature Selector
- [`create_mlp_learner()`](https://kstawiski.github.io/OmicSelector/reference/create_mlp_learner.md)
  : Create MLP Learner via mlr3torch
- [`create_omic_cache()`](https://kstawiski.github.io/OmicSelector/reference/create_omic_cache.md)
  : Create Split-Aware Cache
- [`create_report_data()`](https://kstawiski.github.io/OmicSelector/reference/create_report_data.md)
  : Create Report Data Schema
- [`create_stability_ensemble()`](https://kstawiski.github.io/OmicSelector/reference/create_stability_ensemble.md)
  : Create a Stability Ensemble with Presets
- [`create_ws_pipeop()`](https://kstawiski.github.io/OmicSelector/reference/create_ws_pipeop.md)
  : Apply Within-Sample Normalization to OmicPipeline
- [`cross-platform`](https://kstawiski.github.io/OmicSelector/reference/cross-platform.md)
  : Cross-Platform Transfer Learning Utilities
- [`cross_platform_transfer()`](https://kstawiski.github.io/OmicSelector/reference/cross_platform_transfer.md)
  : Convenience Wrapper for Cross-Platform Transfer
- [`decompose_brier()`](https://kstawiski.github.io/OmicSelector/reference/decompose_brier.md)
  : Decompose Brier Score
- [`encode_batch()`](https://kstawiski.github.io/OmicSelector/reference/encode_batch.md)
  : Batch-Encode an Expression Matrix with Any Encoding
- [`encode_corr_grid()`](https://kstawiski.github.io/OmicSelector/reference/encode_corr_grid.md)
  : Encode with Correlation-Ordered Grid
- [`encode_deepinsight()`](https://kstawiski.github.io/OmicSelector/reference/encode_deepinsight.md)
  : Encode with DeepInsight Layout
- [`encode_ratio_image()`](https://kstawiski.github.io/OmicSelector/reference/encode_ratio_image.md)
  : Create Pairwise Log-Ratio Image
- [`encode_simple_grid()`](https://kstawiski.github.io/OmicSelector/reference/encode_simple_grid.md)
  : Simple Row-Major Grid Encoding
- [`evaluate_feature_selection()`](https://kstawiski.github.io/OmicSelector/reference/evaluate_feature_selection.md)
  : Evaluate Feature Selection for Leakage
- [`export_bundle()`](https://kstawiski.github.io/OmicSelector/reference/export_bundle.md)
  : Create Complete Export Bundle
- [`export_mlr3torch_checkpoint()`](https://kstawiski.github.io/OmicSelector/reference/export_mlr3torch_checkpoint.md)
  : Export mlr3torch Checkpoint
- [`export_omicfit_checkpoint()`](https://kstawiski.github.io/OmicSelector/reference/export_omicfit_checkpoint.md)
  : Export mlr3torch Checkpoint from OmicFit
- [`export_onnx()`](https://kstawiski.github.io/OmicSelector/reference/export_onnx.md)
  : Export Model to ONNX Format
- [`export_vetiver()`](https://kstawiski.github.io/OmicSelector/reference/export_vetiver.md)
  : Export Model as Vetiver for Deployment
- [`extract_features_from_resample()`](https://kstawiski.github.io/OmicSelector/reference/extract_features_from_resample.md)
  : Extract Features from All Folds in ResampleResult
- [`extract_selected_features()`](https://kstawiski.github.io/OmicSelector/reference/extract_selected_features.md)
  : Extract Selected Features from Trained GraphLearner
- [`feature_importance()`](https://kstawiski.github.io/OmicSelector/reference/feature_importance.md)
  : Compute Permutation Feature Importance
- [`finetune_mlr3torch_checkpoint()`](https://kstawiski.github.io/OmicSelector/reference/finetune_mlr3torch_checkpoint.md)
  : Fine-tune from mlr3torch Checkpoint
- [`finetune_omicfit_checkpoint()`](https://kstawiski.github.io/OmicSelector/reference/finetune_omicfit_checkpoint.md)
  : Fine-tune OmicFit from mlr3torch Checkpoint
- [`fit_compositional_mahalanobis()`](https://kstawiski.github.io/OmicSelector/reference/fit_compositional_mahalanobis.md)
  : Fit robust Mahalanobis detector on compositional log-ratio
  coordinates
- [`fit_conformal_anomaly()`](https://kstawiski.github.io/OmicSelector/reference/fit_conformal_anomaly.md)
  : Fit a conformal anomaly detector from a Tier R healthy reference
  cohort
- [`fit_corr_order()`](https://kstawiski.github.io/OmicSelector/reference/fit_corr_order.md)
  : Correlation-Ordered Grid Encoding
- [`fit_deepinsight()`](https://kstawiski.github.io/OmicSelector/reference/fit_deepinsight.md)
  : Fit DeepInsight Feature Layout
- [`fit_dominance_threshold()`](https://kstawiski.github.io/OmicSelector/reference/fit_dominance_threshold.md)
  : Calibrate a cohort-specific dominance threshold
- [`fit_frozen_quantile()`](https://kstawiski.github.io/OmicSelector/reference/fit_frozen_quantile.md)
  : Fit a frozen monotone quantile calibrator from training data
- [`fit_frozen_ruv()`](https://kstawiski.github.io/OmicSelector/reference/fit_frozen_ruv.md)
  : Fit a frozen RUV factor model from training data
- [`fit_hemolysis_prefilter()`](https://kstawiski.github.io/OmicSelector/reference/fit_hemolysis_prefilter.md)
  : Fit a Hemolysis Pre-Filter
- [`fit_hemolysis_rr()`](https://kstawiski.github.io/OmicSelector/reference/fit_hemolysis_rr.md)
  : Fit robust-regression hemolysis nuisance model from training
  controls
- [`fit_isolation_forest_logratio()`](https://kstawiski.github.io/OmicSelector/reference/fit_isolation_forest_logratio.md)
  : Fit a pure-R isolation forest on rCLR log-ratio inputs
- [`fit_isotonic_calibration()`](https://kstawiski.github.io/OmicSelector/reference/fit_isotonic_calibration.md)
  : Isotonic Regression Calibration
- [`fit_logistic_normal_eb()`](https://kstawiski.github.io/OmicSelector/reference/fit_logistic_normal_eb.md)
  : Fit frozen logistic-normal empirical-Bayes prior from training data
- [`fit_platt_scaling()`](https://kstawiski.github.io/OmicSelector/reference/fit_platt_scaling.md)
  : Platt Scaling (Logistic Calibration)
- [`fit_robust_pca_residual()`](https://kstawiski.github.io/OmicSelector/reference/fit_robust_pca_residual.md)
  : Fit a MAD-scaled SVD residual filter for batch correction
- [`fit_temperature_scaling()`](https://kstawiski.github.io/OmicSelector/reference/fit_temperature_scaling.md)
  : Temperature Scaling
- [`frozen-combat`](https://kstawiski.github.io/OmicSelector/reference/frozen-combat.md)
  : Frozen ComBat for Leakage-Free Batch Correction
- [`frozen_combat_correct()`](https://kstawiski.github.io/OmicSelector/reference/frozen_combat_correct.md)
  : Convenience function for frozen ComBat correction
- [`generate_cache_key()`](https://kstawiski.github.io/OmicSelector/reference/generate_cache_key.md)
  : Generate Split-Aware Cache Key
- [`generate_gold_standard()`](https://kstawiski.github.io/OmicSelector/reference/generate_gold_standard.md)
  : Generate Gold Standard Synthetic Dataset
- [`generate_tripod_report()`](https://kstawiski.github.io/OmicSelector/reference/generate_tripod_report.md)
  : Generate TRIPOD+AI Report
- [`get_consensus_features()`](https://kstawiski.github.io/OmicSelector/reference/get_consensus_features.md)
  : Get Consensus Features from Best Signature
- [`get_modality_info()`](https://kstawiski.github.io/OmicSelector/reference/get_modality_info.md)
  : Get Modality Information
- [`get_optimal_params()`](https://kstawiski.github.io/OmicSelector/reference/get_optimal_params.md)
  : Get Optimal Hyperparameters from AutoTuner
- [`get_parallel_status()`](https://kstawiski.github.io/OmicSelector/reference/get_parallel_status.md)
  : Get Current Parallelization Status
- [`get_reliable_shap_features()`](https://kstawiski.github.io/OmicSelector/reference/get_reliable_shap_features.md)
  : Get Reliable SHAP Features
- [`get_selected_features_per_fold()`](https://kstawiski.github.io/OmicSelector/reference/get_selected_features_per_fold.md)
  : Get Selected Features Per Fold
- [`hemolysis-correction`](https://kstawiski.github.io/OmicSelector/reference/hemolysis-correction.md)
  : Hemolysis-Aware Corrections for Biomarker Panels
- [`hemolysis_index_blondal()`](https://kstawiski.github.io/OmicSelector/reference/hemolysis_index_blondal.md)
  : Blondal hemolysis index (log miR-451a - log miR-23a-3p)
- [`hemolysis_proxy_score()`](https://kstawiski.github.io/OmicSelector/reference/hemolysis_proxy_score.md)
  : Hemolysis Proxy Score
- [`image-encodings`](https://kstawiski.github.io/OmicSelector/reference/image-encodings.md)
  : Image Encoding Methods for Biomarker Panels
- [`import_mlr3torch_checkpoint()`](https://kstawiski.github.io/OmicSelector/reference/import_mlr3torch_checkpoint.md)
  : Import mlr3torch Checkpoint
- [`import_omicfit_checkpoint()`](https://kstawiski.github.io/OmicSelector/reference/import_omicfit_checkpoint.md)
  : Import mlr3torch Checkpoint into OmicFit
- [`impute_within_fold()`](https://kstawiski.github.io/OmicSelector/reference/impute_within_fold.md)
  : Within-Fold Median Imputation (Leakage-Free)
- [`interpretability`](https://kstawiski.github.io/OmicSelector/reference/interpretability.md)
  : Model Interpretability for OmicSelector
- [`load_bundle()`](https://kstawiski.github.io/OmicSelector/reference/load_bundle.md)
  : Load Exported Model Bundle
- [`make_autotuner_glmnet()`](https://kstawiski.github.io/OmicSelector/reference/make_autotuner_glmnet.md)
  : Create Bayesian-Optimized AutoTuner for glmnet
- [`make_autotuner_lightgbm()`](https://kstawiski.github.io/OmicSelector/reference/make_autotuner_lightgbm.md)
  : Create Bayesian-Optimized AutoTuner for LightGBM
- [`make_autotuner_ranger()`](https://kstawiski.github.io/OmicSelector/reference/make_autotuner_ranger.md)
  : Create Bayesian-Optimized AutoTuner for Random Forest
- [`make_autotuner_xgboost()`](https://kstawiski.github.io/OmicSelector/reference/make_autotuner_xgboost.md)
  : Create Bayesian-Optimized AutoTuner for XGBoost
- [`make_catboost_learner()`](https://kstawiski.github.io/OmicSelector/reference/make_catboost_learner.md)
  : Create a CatBoost Learner
- [`make_fttransformer_learner()`](https://kstawiski.github.io/OmicSelector/reference/make_fttransformer_learner.md)
  [`make_tabtransformer_learner()`](https://kstawiski.github.io/OmicSelector/reference/make_fttransformer_learner.md)
  : Create an FT-Transformer Learner
- [`make_gnn_learner()`](https://kstawiski.github.io/OmicSelector/reference/make_gnn_learner.md)
  : Create GNN Learner for Pathway-Aware Classification
- [`make_gof_filter()`](https://kstawiski.github.io/OmicSelector/reference/make_gof_filter.md)
  : Create GOF Filter
- [`make_mlp_learner()`](https://kstawiski.github.io/OmicSelector/reference/make_mlp_learner.md)
  : Make MLP Learner (Alias)
- [`make_ratio_image()`](https://kstawiski.github.io/OmicSelector/reference/make_ratio_image.md)
  : Create Pairwise Ratio Image from Expression Vector
- [`make_ratio_images()`](https://kstawiski.github.io/OmicSelector/reference/make_ratio_images.md)
  : Batch-Create Ratio Images from an Expression Matrix
- [`make_tabm_learner()`](https://kstawiski.github.io/OmicSelector/reference/make_tabm_learner.md)
  : Create a TabM Learner
- [`make_tabnet_learner()`](https://kstawiski.github.io/OmicSelector/reference/make_tabnet_learner.md)
  : Create a TabNet Learner
- [`make_tabpfn_learner()`](https://kstawiski.github.io/OmicSelector/reference/make_tabpfn_learner.md)
  : Create a TabPFN Learner
- [`memoize_with_split()`](https://kstawiski.github.io/OmicSelector/reference/memoize_with_split.md)
  : Memoize Function with Split Context
- [`merge_omics_data()`](https://kstawiski.github.io/OmicSelector/reference/merge_omics_data.md)
  : Merge Multi-Omics Data for Analysis
- [`methods-comparison`](https://kstawiski.github.io/OmicSelector/reference/methods-comparison.md)
  : Methods Comparison Utilities
- [`mirna_alias_table()`](https://kstawiski.github.io/OmicSelector/reference/mirna_alias_table.md)
  : Curated miRNA alias lookup table (miRBase v22.1)
- [`model-export`](https://kstawiski.github.io/OmicSelector/reference/model-export.md)
  : Model Export for OmicSelector
- [`multi-omics`](https://kstawiski.github.io/OmicSelector/reference/multi-omics.md)
  : Multi-Omics Support for OmicSelector
- [`noise_augment()`](https://kstawiski.github.io/OmicSelector/reference/noise_augment.md)
  : Gaussian Noise Augmentation
- [`omic_benchmark()`](https://kstawiski.github.io/OmicSelector/reference/omic_benchmark.md)
  : Create benchmark service from OmicPipeline
- [`omic_pipeline()`](https://kstawiski.github.io/OmicSelector/reference/omic_pipeline.md)
  : Quick pipeline creation from data
- [`os_bias_audit()`](https://kstawiski.github.io/OmicSelector/reference/os_bias_audit.md)
  : One-Shot Bias Audit Report
- [`os_bias_audit_report()`](https://kstawiski.github.io/OmicSelector/reference/os_bias_audit_report.md)
  : Plain-Text Bias-Audit Report
- [`os_bias_floor_auc()`](https://kstawiski.github.io/OmicSelector/reference/os_bias_floor_auc.md)
  : Dataset-Identity AUC ("Bias Floor")
- [`os_calibrated_brier()`](https://kstawiski.github.io/OmicSelector/reference/os_calibrated_brier.md)
  : Compute Apparent or Cross-Fitted Calibrated Brier Score
- [`os_claim_gate()`](https://kstawiski.github.io/OmicSelector/reference/os_claim_gate.md)
  : Claim Gate for Panel-Null Benchmarks
- [`os_clustered_bootstrap_auc()`](https://kstawiski.github.io/OmicSelector/reference/os_clustered_bootstrap_auc.md)
  : Clustered Bootstrap CI for AUC
- [`os_conformal_anomaly()`](https://kstawiski.github.io/OmicSelector/reference/os_conformal_anomaly.md)
  : Compute conformal anomaly p-value for new samples
- [`os_conformal_anomaly_score()`](https://kstawiski.github.io/OmicSelector/reference/os_conformal_anomaly_score.md)
  : Conformal Healthy-Reference Anomaly Score
- [`os_covariate_only_auc()`](https://kstawiski.github.io/OmicSelector/reference/os_covariate_only_auc.md)
  : Covariate-Only AUC
- [`os_detect_cross_cohort_duplicates()`](https://kstawiski.github.io/OmicSelector/reference/os_detect_cross_cohort_duplicates.md)
  : Specimen-Duplication Detection Across Cohorts
- [`os_grouped_resample_auc()`](https://kstawiski.github.io/OmicSelector/reference/os_grouped_resample_auc.md)
  : Estimate AUC over Grouped Resampling Folds
- [`os_identifiability_gate()`](https://kstawiski.github.io/OmicSelector/reference/os_identifiability_gate.md)
  : Apply a Fail-Closed Identifiability Gate
- [`os_ktsp_fit()`](https://kstawiski.github.io/OmicSelector/reference/os_ktsp_fit.md)
  : Fit a Top-k Oriented-Pair Panel Scorer
- [`os_log_transform_adaptive()`](https://kstawiski.github.io/OmicSelector/reference/os_log_transform_adaptive.md)
  : Adaptive Pseudocount log2 Transform
- [`os_mahalanobis_score()`](https://kstawiski.github.io/OmicSelector/reference/os_mahalanobis_score.md)
  : Mahalanobis Anomaly Score from a Reference Set
- [`os_make_grouped_stratified_folds()`](https://kstawiski.github.io/OmicSelector/reference/os_make_grouped_stratified_folds.md)
  : Build Grouped Stratified Folds
- [`os_mde_margin_record()`](https://kstawiski.github.io/OmicSelector/reference/os_mde_margin_record.md)
  : Create an MDE/Margin Record
- [`os_null_qc()`](https://kstawiski.github.io/OmicSelector/reference/os_null_qc.md)
  : Validate a Random-Panel Null Benchmark
- [`os_oof_pipeline_compare()`](https://kstawiski.github.io/OmicSelector/reference/os_oof_pipeline_compare.md)
  : Out-of-Fold Pipeline Comparison
- [`os_operating_point_gate()`](https://kstawiski.github.io/OmicSelector/reference/os_operating_point_gate.md)
  : Fail-Closed Operating-Point Gate
- [`os_operating_points()`](https://kstawiski.github.io/OmicSelector/reference/os_operating_points.md)
  : Compute Binary Clinical Operating Points
- [`os_paired_delong()`](https://kstawiski.github.io/OmicSelector/reference/os_paired_delong.md)
  : Paired DeLong Test for AUC Comparison
- [`os_panel_null_benchmark()`](https://kstawiski.github.io/OmicSelector/reference/os_panel_null_benchmark.md)
  : Matched Random-Panel Null Benchmark
- [`os_per_feature_batch_signal()`](https://kstawiski.github.io/OmicSelector/reference/os_per_feature_batch_signal.md)
  : Per-Feature Batch-vs-Case Signal Partitioning
- [`os_plate_median_frozen`](https://kstawiski.github.io/OmicSelector/reference/os_plate_median_frozen.md)
  : Frozen Plate-Median Correction
- [`os_provenance_floor_suite()`](https://kstawiski.github.io/OmicSelector/reference/os_provenance_floor_suite.md)
  : Estimate a Provenance-Only Prediction Floor
- [`os_provenance_preflight()`](https://kstawiski.github.io/OmicSelector/reference/os_provenance_preflight.md)
  : Paper 3 specimen-overlap provenance pre-flight gate
- [`os_required_margin()`](https://kstawiski.github.io/OmicSelector/reference/os_required_margin.md)
  : Required Margin for Matched-Null Claims
- [`os_singscore()`](https://kstawiski.github.io/OmicSelector/reference/os_singscore.md)
  : Score a Direction-Split Panel by Within-Sample Ranks
- [`os_terminal_gate_ledger()`](https://kstawiski.github.io/OmicSelector/reference/os_terminal_gate_ledger.md)
  : Build a Terminal Gate Ledger
- [`os_validate_folds()`](https://kstawiski.github.io/OmicSelector/reference/os_validate_folds.md)
  : Validate Grouped Folds
- [`os_within_provenance_blocks()`](https://kstawiski.github.io/OmicSelector/reference/os_within_provenance_blocks.md)
  : Summarize Within-Provenance Case-Control Blocks
- [`paper3-additional-within-sample`](https://kstawiski.github.io/OmicSelector/reference/paper3-additional-within-sample.md)
  : Paper 3 additional within-sample methods (Module A, P2)
- [`paper3-batch-correction`](https://kstawiski.github.io/OmicSelector/reference/paper3-batch-correction.md)
  : Paper 3 batch-correction methods (Module C)
- [`paper3-hemolysis`](https://kstawiski.github.io/OmicSelector/reference/paper3-hemolysis.md)
  : Paper 3 robust-regression hemolysis correction (Module B)
- [`paper3-matched-null`](https://kstawiski.github.io/OmicSelector/reference/paper3-matched-null.md)
  : Paper 3 matched-null benchmark for panel-vs-random-panel AUC
  inference
- [`paper3-nondetects`](https://kstawiski.github.io/OmicSelector/reference/paper3-nondetects.md)
  : Paper 3 qPCR non-detect imputation (Module B add-on)
- [`paper3-outlier-detection`](https://kstawiski.github.io/OmicSelector/reference/paper3-outlier-detection.md)
  : Paper 3 outlier detection and conformal claim-gating (Module D)
- [`paper3-preprocessing`](https://kstawiski.github.io/OmicSelector/reference/paper3-preprocessing.md)
  : Paper 3 preprocessing utilities
- [`paper3-within-sample`](https://kstawiski.github.io/OmicSelector/reference/paper3-within-sample.md)
  : Paper 3 within-sample compositional methods (Module A, P1)
- [`paper3_bh_fdr_correct_blocked()`](https://kstawiski.github.io/OmicSelector/reference/paper3_bh_fdr_correct_blocked.md)
  : Block-aware two-stage BH-FDR for specimen-shared cohort clusters
- [`paper3_bh_fdr_correct_matched_null()`](https://kstawiski.github.io/OmicSelector/reference/paper3_bh_fdr_correct_matched_null.md)
  : BH-FDR correction within a matched-null modality family
- [`paper3_hanley_mcneil_auc_ci()`](https://kstawiski.github.io/OmicSelector/reference/paper3_hanley_mcneil_auc_ci.md)
  : Hanley-McNeil 1982 AUC confidence interval
- [`paper3_holm_correct_familywise()`](https://kstawiski.github.io/OmicSelector/reference/paper3_holm_correct_familywise.md)
  : Holm correction for family-wise method contrasts
- [`paper3_matched_null_benchmark()`](https://kstawiski.github.io/OmicSelector/reference/paper3_matched_null_benchmark.md)
  : Paper 3 matched-null benchmark for within-sample miRNA panel scoring
- [`paper3_matched_null_benchmark_cv()`](https://kstawiski.github.io/OmicSelector/reference/paper3_matched_null_benchmark_cv.md)
  : 5-fold nested-CV matched-null benchmark
- [`parallel`](https://kstawiski.github.io/OmicSelector/reference/parallel.md)
  : Parallelization Support for OmicSelector
- [`partial_dependence()`](https://kstawiski.github.io/OmicSelector/reference/partial_dependence.md)
  : Compute Partial Dependence
- [`plot_signature_tradeoffs()`](https://kstawiski.github.io/OmicSelector/reference/plot_signature_tradeoffs.md)
  : Plot Signature Selection Trade-offs
- [`plot_xai_importance()`](https://kstawiski.github.io/OmicSelector/reference/plot_xai_importance.md)
  : Plot XAI Feature Importance
- [`predict(`*`<os_ktsp_model>`*`)`](https://kstawiski.github.io/OmicSelector/reference/predict.os_ktsp_model.md)
  : Predict from a Top-k Oriented-Pair Panel Scorer
- [`predict_clr_mlp()`](https://kstawiski.github.io/OmicSelector/reference/predict_clr_mlp.md)
  : Predict with a CLR + MLP Classifier
- [`predict_codacore()`](https://kstawiski.github.io/OmicSelector/reference/predict_codacore.md)
  : Predict with a CoDaCoRe Classifier
- [`predict_weighted_clr_mlp()`](https://kstawiski.github.io/OmicSelector/reference/predict_weighted_clr_mlp.md)
  : Predict with a Weighted CLR + MLP Classifier
- [`preprocess_inverse_log()`](https://kstawiski.github.io/OmicSelector/reference/preprocess_inverse_log.md)
  : Inverse-log preprocessor for pre-log-transformed microarray deposits
- [`print(`*`<CalibrationResult>`*`)`](https://kstawiski.github.io/OmicSelector/reference/print.CalibrationResult.md)
  : Print method for CalibrationResult
- [`print(`*`<CorrelationCheck>`*`)`](https://kstawiski.github.io/OmicSelector/reference/print.CorrelationCheck.md)
  : Print Correlation Check
- [`print(`*`<FeatureImportance>`*`)`](https://kstawiski.github.io/OmicSelector/reference/print.FeatureImportance.md)
  : Print Feature Importance
- [`print(`*`<NestedCVResult>`*`)`](https://kstawiski.github.io/OmicSelector/reference/print.NestedCVResult.md)
  : Print method for NestedCVResult
- [`print(`*`<NogueiraStability>`*`)`](https://kstawiski.github.io/OmicSelector/reference/print.NogueiraStability.md)
  : Print method for NogueiraStability
- [`print(`*`<OmicBenchmarkResult>`*`)`](https://kstawiski.github.io/OmicSelector/reference/print.OmicBenchmarkResult.md)
  : Print method for OmicBenchmarkResult
- [`print(`*`<OmicFit>`*`)`](https://kstawiski.github.io/OmicSelector/reference/print.OmicFit.md)
  : Print method for OmicFit
- [`print(`*`<OmicsInput>`*`)`](https://kstawiski.github.io/OmicSelector/reference/print.OmicsInput.md)
  : Print method for OmicsInput
- [`print(`*`<ReportData>`*`)`](https://kstawiski.github.io/OmicSelector/reference/print.ReportData.md)
  : Print method for ReportData
- [`print(`*`<SignatureSelectionResult>`*`)`](https://kstawiski.github.io/OmicSelector/reference/print.SignatureSelectionResult.md)
  : Print Method for Signature Selection Result
- [`print_shap_warnings()`](https://kstawiski.github.io/OmicSelector/reference/print_shap_warnings.md)
  : Print SHAP Warnings Report
- [`print_xai_summary()`](https://kstawiski.github.io/OmicSelector/reference/print_xai_summary.md)
  : Print XAI Summary
- [`qpcr_nondetect_impute()`](https://kstawiski.github.io/OmicSelector/reference/qpcr_nondetect_impute.md)
  : Bayesian hierarchical imputation of qPCR non-detects
- [`qpcr_nondetect_lod_fallback()`](https://kstawiski.github.io/OmicSelector/reference/qpcr_nondetect_lod_fallback.md)
  : Limit-of-detection fallback imputation for qPCR non-detects
- [`ratio-image-cnn`](https://kstawiski.github.io/OmicSelector/reference/ratio-image-cnn.md)
  : Ratio Image CNN for Biomarker Panel Classification
- [`register_coda_feature_selection_filters()`](https://kstawiski.github.io/OmicSelector/reference/register_coda_feature_selection_filters.md)
  : Register CoDA Feature Selection Filters
- [`register_gof_filters()`](https://kstawiski.github.io/OmicSelector/reference/register_gof_filters.md)
  : Register GOF Filters in mlr3
- [`reliability_diagram_data()`](https://kstawiski.github.io/OmicSelector/reference/reliability_diagram_data.md)
  : Create Reliability Diagram Data
- [`report`](https://kstawiski.github.io/OmicSelector/reference/report.md)
  : TRIPOD+AI Report Generation for OmicSelector
- [`reset_parallel()`](https://kstawiski.github.io/OmicSelector/reference/reset_parallel.md)
  : Reset Parallelization to Sequential
- [`resolve_mirna_aliases()`](https://kstawiski.github.io/OmicSelector/reference/resolve_mirna_aliases.md)
  : Resolve miRNA identifiers to a canonical namespace
- [`run_bayesian_benchmark()`](https://kstawiski.github.io/OmicSelector/reference/run_bayesian_benchmark.md)
  : Run Bayesian Optimization Benchmark
- [`run_dl_benchmark()`](https://kstawiski.github.io/OmicSelector/reference/run_dl_benchmark.md)
  : Deep Learning Benchmark
- [`safe-evaluation`](https://kstawiski.github.io/OmicSelector/reference/safe-evaluation.md)
  : Safe Evaluation Utilities for Biomarker Validation
- [`safe-preprocessing`](https://kstawiski.github.io/OmicSelector/reference/safe-preprocessing.md)
  : Safe Preprocessing Utilities for Cross-Validation
- [`safe_auc()`](https://kstawiski.github.io/OmicSelector/reference/safe_auc.md)
  : Safe AUC with Confidence Interval
- [`safe_roc()`](https://kstawiski.github.io/OmicSelector/reference/safe_roc.md)
  : Safe ROC Computation (Direction-Guarded)
- [`select_best_signature()`](https://kstawiski.github.io/OmicSelector/reference/select_best_signature.md)
  : Select Best Biomarker Signature from Nested CV Results
- [`setup_parallel()`](https://kstawiski.github.io/OmicSelector/reference/setup_parallel.md)
  : Configure Parallelization for OmicSelector
- [`shap_values()`](https://kstawiski.github.io/OmicSelector/reference/shap_values.md)
  : Compute SHAP-like Values
- [`shap_warnings`](https://kstawiski.github.io/OmicSelector/reference/shap_warnings.md)
  : Correlation-Aware SHAP Interpretation
- [`signature-selection`](https://kstawiski.github.io/OmicSelector/reference/signature-selection.md)
  : Signature Selection: Multi-Objective Best Biomarker Selection
- [`smote_augment()`](https://kstawiski.github.io/OmicSelector/reference/smote_augment.md)
  : SMOTE Augmentation for Omics Data
- [`stability`](https://kstawiski.github.io/OmicSelector/reference/stability.md)
  : Nogueira Stability Index for Feature Selection
- [`stack_omics()`](https://kstawiski.github.io/OmicSelector/reference/stack_omics.md)
  : Create Multi-Omics Stacked Ensemble (Convenience Function)
- [`standardize_within_fold()`](https://kstawiski.github.io/OmicSelector/reference/standardize_within_fold.md)
  : Within-Fold Standardization (Leakage-Free)
- [`tabddpm_generate()`](https://kstawiski.github.io/OmicSelector/reference/tabddpm_generate.md)
  : TabDDPM Synthetic Data Generator
- [`test_noninferiority()`](https://kstawiski.github.io/OmicSelector/reference/test_noninferiority.md)
  : Paired DeLong Non-Inferiority Test
- [`torch-checkpoint`](https://kstawiski.github.io/OmicSelector/reference/torch-checkpoint.md)
  : mlr3torch Checkpoint Utilities
- [`torch-learners`](https://kstawiski.github.io/OmicSelector/reference/torch-learners.md)
  : mlr3torch Learner Integration for OmicSelector
- [`train_clr_mlp()`](https://kstawiski.github.io/OmicSelector/reference/train_clr_mlp.md)
  : Train a CLR + MLP Classifier
- [`train_codacore()`](https://kstawiski.github.io/OmicSelector/reference/train_codacore.md)
  : Train a CoDaCoRe Classifier
- [`train_ratio_cnn()`](https://kstawiski.github.io/OmicSelector/reference/train_ratio_cnn.md)
  : Train Ratio Image CNN
- [`train_ratio_cnn_multiseed()`](https://kstawiski.github.io/OmicSelector/reference/train_ratio_cnn_multiseed.md)
  : Train Ratio CNN Across Three Default Seeds
- [`train_weighted_clr_mlp()`](https://kstawiski.github.io/OmicSelector/reference/train_weighted_clr_mlp.md)
  : Train a Weighted CLR + MLP Classifier
- [`validate_omics_input()`](https://kstawiski.github.io/OmicSelector/reference/validate_omics_input.md)
  : Validate Multi-Omics Input
- [`validate_synthetic()`](https://kstawiski.github.io/OmicSelector/reference/validate_synthetic.md)
  : Validate Synthetic Data Quality
- [`weighted_clr_transform()`](https://kstawiski.github.io/OmicSelector/reference/weighted_clr_transform.md)
  : Weighted Centered Log-Ratio Transformation
- [`with_parallel()`](https://kstawiski.github.io/OmicSelector/reference/with_parallel.md)
  : With Parallel Scope
- [`within-sample`](https://kstawiski.github.io/OmicSelector/reference/within-sample.md)
  : Within-Sample Normalization for Biomarker Panels
- [`ws_alr_pivot()`](https://kstawiski.github.io/OmicSelector/reference/ws_alr_pivot.md)
  : Additive log-ratio with a frozen pivot pool
- [`ws_balance_ilr()`](https://kstawiski.github.io/OmicSelector/reference/ws_balance_ilr.md)
  : Isometric log-ratio balances on a frozen partition tree
- [`ws_default_pivot_pool()`](https://kstawiski.github.io/OmicSelector/reference/ws_default_pivot_pool.md)
  : Default ALR pivot pool for circulating-miRNA panels (v1)
- [`ws_default_sbp()`](https://kstawiski.github.io/OmicSelector/reference/ws_default_sbp.md)
  : Default circulating-miRNA sequential binary partition (v1)
- [`ws_dominance_flag()`](https://kstawiski.github.io/OmicSelector/reference/ws_dominance_flag.md)
  : Flag samples whose dominance score exceeds a calibrated threshold
- [`ws_dominance_score()`](https://kstawiski.github.io/OmicSelector/reference/ws_dominance_score.md)
  : Panel-internal dominance score (within-sample QC)
- [`ws_logratio()`](https://kstawiski.github.io/OmicSelector/reference/ws_logratio.md)
  : Within-Sample Pairwise Log-Ratio Features
- [`ws_mad_logratio()`](https://kstawiski.github.io/OmicSelector/reference/ws_mad_logratio.md)
  : Median-centred log-ratio with optional MAD scaling
- [`ws_minmax()`](https://kstawiski.github.io/OmicSelector/reference/ws_minmax.md)
  : Within-Sample Min-Max Normalization
- [`ws_perturbation_benchmark()`](https://kstawiski.github.io/OmicSelector/reference/ws_perturbation_benchmark.md)
  : Within-Sample Perturbation Benchmark
- [`ws_rank()`](https://kstawiski.github.io/OmicSelector/reference/ws_rank.md)
  : Within-Sample Rank Normalization
- [`ws_ratio_image()`](https://kstawiski.github.io/OmicSelector/reference/ws_ratio_image.md)
  : Within-Sample Pairwise Ratio Image
- [`ws_rclr_trimmed()`](https://kstawiski.github.io/OmicSelector/reference/ws_rclr_trimmed.md)
  : Robust trimmed centred log-ratio for compositional miRNA panels
- [`ws_zscore()`](https://kstawiski.github.io/OmicSelector/reference/ws_zscore.md)
  : Within-Sample Z-Score Normalization
- [`xai_correlations()`](https://kstawiski.github.io/OmicSelector/reference/xai_correlations.md)
  : Compute Correlation Diagnostics for Features
- [`xai_explainer_mlr3()`](https://kstawiski.github.io/OmicSelector/reference/xai_explainer_mlr3.md)
  : Create DALEX Explainer from mlr3 Learner
- [`xai_importance()`](https://kstawiski.github.io/OmicSelector/reference/xai_importance.md)
  : Compute Permutation Feature Importance
- [`xai_pdp()`](https://kstawiski.github.io/OmicSelector/reference/xai_pdp.md)
  : Compute Partial Dependence Plots
- [`xai_pipeline()`](https://kstawiski.github.io/OmicSelector/reference/xai_pipeline.md)
  : Run Complete XAI Pipeline
- [`xai_shap()`](https://kstawiski.github.io/OmicSelector/reference/xai_shap.md)
  : Compute SHAP Values for Observations
