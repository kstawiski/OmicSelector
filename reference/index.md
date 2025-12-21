# Package index

## Core Pipeline

Main classes for building zero-leakage ML pipelines

- [`OmicPipeline`](https://biostat.umed.pl/OmicSelector/reference/OmicPipeline.md)
  : OmicPipeline: Zero-Leakage Feature Selection Pipeline
- [`omic_pipeline()`](https://biostat.umed.pl/OmicSelector/reference/omic_pipeline.md)
  : Quick pipeline creation from data
- [`BenchmarkService`](https://biostat.umed.pl/OmicSelector/reference/BenchmarkService.md)
  : BenchmarkService: Nested Cross-Validation with Zero Leakage
- [`omic_benchmark()`](https://biostat.umed.pl/OmicSelector/reference/omic_benchmark.md)
  : Create benchmark service from OmicPipeline

## Signature Selection

Multi-objective biomarker signature selection

- [`select_best_signature()`](https://biostat.umed.pl/OmicSelector/reference/select_best_signature.md)
  : Select Best Biomarker Signature from Nested CV Results
- [`get_consensus_features()`](https://biostat.umed.pl/OmicSelector/reference/get_consensus_features.md)
  : Get Consensus Features from Best Signature
- [`get_selected_features_per_fold()`](https://biostat.umed.pl/OmicSelector/reference/get_selected_features_per_fold.md)
  : Get Selected Features Per Fold

## Stability Analysis

Feature selection stability metrics

- [`compute_nogueira_stability()`](https://biostat.umed.pl/OmicSelector/reference/compute_nogueira_stability.md)
  : Compute Nogueira Stability Index
- [`compute_stability_from_resample()`](https://biostat.umed.pl/OmicSelector/reference/compute_stability_from_resample.md)
  : Compute Stability from ResampleResult
- [`extract_features_from_resample()`](https://biostat.umed.pl/OmicSelector/reference/extract_features_from_resample.md)
  : Extract Features from All Folds in ResampleResult
- [`extract_selected_features()`](https://biostat.umed.pl/OmicSelector/reference/extract_selected_features.md)
  : Extract Selected Features from Trained GraphLearner

## Batch Correction

FrozenComBat for proper batch effect handling

- [`FrozenComBat`](https://biostat.umed.pl/OmicSelector/reference/FrozenComBat.md)
  : FrozenComBat R6 Class
- [`frozen_combat_correct()`](https://biostat.umed.pl/OmicSelector/reference/frozen_combat_correct.md)
  : Convenience function for frozen ComBat correction
- [`create_frozen_combat_pipeop()`](https://biostat.umed.pl/OmicSelector/reference/create_frozen_combat_pipeop.md)
  : Create a Frozen ComBat PipeOp for mlr3pipelines

## Calibration

Probability calibration and diagnostics

- [`fit_platt_scaling()`](https://biostat.umed.pl/OmicSelector/reference/fit_platt_scaling.md)
  : Platt Scaling (Logistic Calibration)
- [`fit_isotonic_calibration()`](https://biostat.umed.pl/OmicSelector/reference/fit_isotonic_calibration.md)
  : Isotonic Regression Calibration
- [`fit_temperature_scaling()`](https://biostat.umed.pl/OmicSelector/reference/fit_temperature_scaling.md)
  : Temperature Scaling
- [`compute_ece()`](https://biostat.umed.pl/OmicSelector/reference/compute_ece.md)
  : Compute Expected Calibration Error (ECE)
- [`decompose_brier()`](https://biostat.umed.pl/OmicSelector/reference/decompose_brier.md)
  : Decompose Brier Score
- [`calibration_summary()`](https://biostat.umed.pl/OmicSelector/reference/calibration_summary.md)
  : Calibration Summary for Model Results
- [`reliability_diagram_data()`](https://biostat.umed.pl/OmicSelector/reference/reliability_diagram_data.md)
  : Create Reliability Diagram Data

## Interpretability

Model interpretation and feature importance

- [`create_explainer()`](https://biostat.umed.pl/OmicSelector/reference/create_explainer.md)
  : Create Model Explainer
- [`shap_values()`](https://biostat.umed.pl/OmicSelector/reference/shap_values.md)
  : Compute SHAP-like Values
- [`feature_importance()`](https://biostat.umed.pl/OmicSelector/reference/feature_importance.md)
  : Compute Permutation Feature Importance
- [`partial_dependence()`](https://biostat.umed.pl/OmicSelector/reference/partial_dependence.md)
  : Compute Partial Dependence
- [`check_feature_correlations()`](https://biostat.umed.pl/OmicSelector/reference/check_feature_correlations.md)
  : Check Feature Correlations

## Phase 5: GOF Filters

Goodness-of-fit filters for sparse/zero-inflated omics data

- [`FilterGOF_KS`](https://biostat.umed.pl/OmicSelector/reference/FilterGOF_KS.md)
  : Kolmogorov-Smirnov GOF Filter
- [`FilterHurdle`](https://biostat.umed.pl/OmicSelector/reference/FilterHurdle.md)
  : Hurdle Filter for Zero-Inflated Data
- [`FilterZeroProp`](https://biostat.umed.pl/OmicSelector/reference/FilterZeroProp.md)
  : Zero-Proportion Filter
- [`make_gof_filter()`](https://biostat.umed.pl/OmicSelector/reference/make_gof_filter.md)
  : Create GOF Filter
- [`compare_gof_filters()`](https://biostat.umed.pl/OmicSelector/reference/compare_gof_filters.md)
  : Compare GOF Filters on Task
- [`register_gof_filters()`](https://biostat.umed.pl/OmicSelector/reference/register_gof_filters.md)
  : Register GOF Filters in mlr3

## Phase 5: Bayesian Tuning

Bayesian hyperparameter optimization with mlr3mbo

- [`make_autotuner_glmnet()`](https://biostat.umed.pl/OmicSelector/reference/make_autotuner_glmnet.md)
  : Create Bayesian-Optimized AutoTuner for glmnet
- [`make_autotuner_xgboost()`](https://biostat.umed.pl/OmicSelector/reference/make_autotuner_xgboost.md)
  : Create Bayesian-Optimized AutoTuner for XGBoost
- [`make_autotuner_ranger()`](https://biostat.umed.pl/OmicSelector/reference/make_autotuner_ranger.md)
  : Create Bayesian-Optimized AutoTuner for Random Forest
- [`make_autotuner_lightgbm()`](https://biostat.umed.pl/OmicSelector/reference/make_autotuner_lightgbm.md)
  : Create Bayesian-Optimized AutoTuner for LightGBM
- [`get_optimal_params()`](https://biostat.umed.pl/OmicSelector/reference/get_optimal_params.md)
  : Get Optimal Hyperparameters from AutoTuner
- [`run_bayesian_benchmark()`](https://biostat.umed.pl/OmicSelector/reference/run_bayesian_benchmark.md)
  : Run Bayesian Optimization Benchmark

## Phase 5: AutoXAI

DALEX-based interpretability with correlation warnings

- [`xai_pipeline()`](https://biostat.umed.pl/OmicSelector/reference/xai_pipeline.md)
  : Run Complete XAI Pipeline
- [`xai_explainer_mlr3()`](https://biostat.umed.pl/OmicSelector/reference/xai_explainer_mlr3.md)
  : Create DALEX Explainer from mlr3 Learner
- [`xai_importance()`](https://biostat.umed.pl/OmicSelector/reference/xai_importance.md)
  : Compute Permutation Feature Importance
- [`xai_correlations()`](https://biostat.umed.pl/OmicSelector/reference/xai_correlations.md)
  : Compute Correlation Diagnostics for Features
- [`xai_pdp()`](https://biostat.umed.pl/OmicSelector/reference/xai_pdp.md)
  : Compute Partial Dependence Plots
- [`xai_shap()`](https://biostat.umed.pl/OmicSelector/reference/xai_shap.md)
  : Compute SHAP Values for Observations
- [`plot_xai_importance()`](https://biostat.umed.pl/OmicSelector/reference/plot_xai_importance.md)
  : Plot XAI Feature Importance
- [`print_xai_summary()`](https://biostat.umed.pl/OmicSelector/reference/print_xai_summary.md)
  : Print XAI Summary

## Phase 5: Stability Ensemble

Bootstrap-based feature selection stability

- [`create_stability_ensemble()`](https://biostat.umed.pl/OmicSelector/reference/create_stability_ensemble.md)
  : Create a Stability Ensemble with Presets
- [`StabilityEnsemble`](https://biostat.umed.pl/OmicSelector/reference/StabilityEnsemble.md)
  : Stability-Based Ensemble Selection
- [`SequentialSelector`](https://biostat.umed.pl/OmicSelector/reference/SequentialSelector.md)
  : Hybrid Sequential Feature Selection (HSFS)

## Phase 5: Synthetic Data

Data augmentation and synthetic data generation

- [`smote_augment()`](https://biostat.umed.pl/OmicSelector/reference/smote_augment.md)
  : SMOTE Augmentation for Omics Data
- [`noise_augment()`](https://biostat.umed.pl/OmicSelector/reference/noise_augment.md)
  : Gaussian Noise Augmentation
- [`validate_synthetic()`](https://biostat.umed.pl/OmicSelector/reference/validate_synthetic.md)
  : Validate Synthetic Data Quality
- [`balance_classes()`](https://biostat.umed.pl/OmicSelector/reference/balance_classes.md)
  : Create Balanced Training Set
- [`tabddpm_generate()`](https://biostat.umed.pl/OmicSelector/reference/tabddpm_generate.md)
  : TabDDPM Synthetic Data Generator

## Phase 5: Deep Learning

Deep learning models for omics (requires torch)

- [`create_mlp_learner()`](https://biostat.umed.pl/OmicSelector/reference/create_mlp_learner.md)
  : Create MLP Learner via mlr3torch
- [`make_mlp_learner()`](https://biostat.umed.pl/OmicSelector/reference/make_mlp_learner.md)
  : Create Omics-Optimized MLP Learner
- [`make_tabtransformer_learner()`](https://biostat.umed.pl/OmicSelector/reference/make_tabtransformer_learner.md)
  : Create TabTransformer Learner
- [`make_gnn_learner()`](https://biostat.umed.pl/OmicSelector/reference/make_gnn_learner.md)
  : Create GNN Learner for Pathway-Aware Classification
- [`build_correlation_adjacency()`](https://biostat.umed.pl/OmicSelector/reference/build_correlation_adjacency.md)
  : Create Correlation-Based Adjacency for GNN
- [`run_dl_benchmark()`](https://biostat.umed.pl/OmicSelector/reference/run_dl_benchmark.md)
  : Deep Learning Benchmark
- [`check_dl_availability()`](https://biostat.umed.pl/OmicSelector/reference/check_dl_availability.md)
  : Check Deep Learning Availability

## Multi-Omics

Multi-omics late integration

- [`merge_omics_data()`](https://biostat.umed.pl/OmicSelector/reference/merge_omics_data.md)
  : Merge Multi-Omics Data for Analysis
- [`get_modality_info()`](https://biostat.umed.pl/OmicSelector/reference/get_modality_info.md)
  : Get Modality Information
- [`validate_omics_input()`](https://biostat.umed.pl/OmicSelector/reference/validate_omics_input.md)
  : Validate Multi-Omics Input
- [`stack_omics()`](https://biostat.umed.pl/OmicSelector/reference/stack_omics.md)
  : Create Multi-Omics Stacked Ensemble (Convenience Function)

## Model Export & Reporting

Model deployment and TRIPOD reporting

- [`generate_tripod_report()`](https://biostat.umed.pl/OmicSelector/reference/generate_tripod_report.md)
  : Generate TRIPOD+AI Report
- [`create_report_data()`](https://biostat.umed.pl/OmicSelector/reference/create_report_data.md)
  : Create Report Data Schema

## Utilities

Helper functions for parallel processing and caching

- [`setup_parallel()`](https://biostat.umed.pl/OmicSelector/reference/setup_parallel.md)
  : Configure Parallelization for OmicSelector
- [`with_parallel()`](https://biostat.umed.pl/OmicSelector/reference/with_parallel.md)
  : With Parallel Scope
- [`get_parallel_status()`](https://biostat.umed.pl/OmicSelector/reference/get_parallel_status.md)
  : Get Current Parallelization Status
- [`reset_parallel()`](https://biostat.umed.pl/OmicSelector/reference/reset_parallel.md)
  : Reset Parallelization to Sequential
- [`create_omic_cache()`](https://biostat.umed.pl/OmicSelector/reference/create_omic_cache.md)
  : Create Split-Aware Cache
- [`cached_filter()`](https://biostat.umed.pl/OmicSelector/reference/cached_filter.md)
  : Cached Filter Computation
- [`cache_stats()`](https://biostat.umed.pl/OmicSelector/reference/cache_stats.md)
  : Get Cache Statistics
- [`clear_cache()`](https://biostat.umed.pl/OmicSelector/reference/clear_cache.md)
  : Clear OmicSelector Cache

## Legacy Functions

Functions from OmicSelector 1.x (use with caution)

- [`list_deprecated_functions()`](https://biostat.umed.pl/OmicSelector/reference/list_deprecated_functions.md)
  : List all deprecated functions
- [`OmicSelector_OmicSelector_wrapper()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_OmicSelector_wrapper.md)
  : OmicSelector_OmicSelector (DEPRECATED)
- [`OmicSelector_PCA()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_PCA.md)
  : OmicSelector_PCA
- [`OmicSelector_PCA_3D()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_PCA_3D.md)
  : OmicSelector_PCA_3D
- [`OmicSelector_benchmark()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_benchmark.md)
  : OmicSelector_benchmark
- [`OmicSelector_best_signature_de()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_best_signature_de.md)
  : OmicSelector_best_signature_de
- [`OmicSelector_best_signature_proposals()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_best_signature_proposals.md)
  : OmicSelector_best_signature_proposals
- [`OmicSelector_best_signature_proposals_meta11()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_best_signature_proposals_meta11.md)
  : OmicSelector_best_signature_proposals_meta11
- [`OmicSelector_combat()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_combat.md)
  : OmicSelector_combat
- [`OmicSelector_correct_miRNA_names()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_correct_miRNA_names.md)
  : OmicSelector_correct_miRNA_names
- [`OmicSelector_correlation_plot()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_correlation_plot.md)
  : OmicSelector_correlation_plot
- [`OmicSelector_counts_to_log10tpm()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_counts_to_log10tpm.md)
  : OmicSelector_counts_to_log10tpm
- [`OmicSelector_create_formula()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_create_formula.md)
  : OmicSelector_create_formula
- [`OmicSelector_differential_expression_ttest()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_differential_expression_ttest.md)
  : OmicSelector_differential_expression_ttest
- [`OmicSelector_diverge_color()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_diverge_color.md)
  : OmicSelector_diverge_color
- [`OmicSelector_docker.update_progress()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_docker.update_progress.md)
  : OmicSelector_docker.update_progress
- [`OmicSelector_download_tissue_miRNA_data_from_TCGA()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_download_tissue_miRNA_data_from_TCGA.md)
  : OmicSelector_download_tissue_miRNA_data_from_TCGA
- [`OmicSelector_get_benchmark()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_get_benchmark.md)
  : OmicSelector_get_benchmark
- [`OmicSelector_get_benchmark_methods()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_get_benchmark_methods.md)
  : OmicSelector_get_benchmark_methods
- [`OmicSelector_get_features_from_benchmark()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_get_features_from_benchmark.md)
  : OmicSelector_get_features_from_benchmark
- [`OmicSelector_heatmap.3()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_heatmap.3.md)
  : OmicSelector_heatmap.3
- [`OmicSelector_heatmap()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_heatmap.md)
  : OmicSelector_heatmap
- [`OmicSelector_iteratedRFE()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_iteratedRFE.md)
  : OmicSelector_iteratedRFE
- [`OmicSelector_iteratedRFE_deprecated()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_iteratedRFE_deprecated.md)
  : Wrapper for deprecated OmicSelector_iteratedRFE
- [`OmicSelector_load_datamix()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_load_datamix.md)
  : OmicSelector_load_datamix
- [`OmicSelector_load_extension()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_load_extension.md)
  : OmicSelector_load_extension
- [`OmicSelector_log()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_log.md)
  : OmicSelector_log
- [`OmicSelector_merge_formulas()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_merge_formulas.md)
  : OmicSelector_merge_formulas
- [`OmicSelector_myclust()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_myclust.md)
  : OmicSelector_myclust
- [`OmicSelector_mydist()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_mydist.md)
  : OmicSelector_mydist
- [`OmicSelector_prepare_split()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_prepare_split.md)
  : OmicSelector_prepare_split
- [`OmicSelector_process_tissue_miRNA_TCGA()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_process_tissue_miRNA_TCGA.md)
  : OmicSelector_process_tissue_miRNA_TCGA
- [`OmicSelector_profileplot()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_profileplot.md)
  : OmicSelector_profileplot
- [`OmicSelector_propensity_score_matching()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_propensity_score_matching.md)
  : OmicSelector_propensity_score_matching
- [`OmicSelector_setup()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_setup.md)
  : OmicSelector_setup
- [`OmicSelector_signature_overlap()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_signature_overlap.md)
  : OmicSelector_signature_overlap
- [`OmicSelector_table()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_table.md)
  : OmicSelector_table
- [`OmicSelector_vulcano_plot()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_vulcano_plot.md)
  : OmicSelector_vulcano_plot
- [`OmicSelector_xgboost()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_xgboost.md)
  : OmicSelector_xgboost

## Internal Functions

Internal functions and helpers (exported for advanced use)

- [`AutoXAI`](https://biostat.umed.pl/OmicSelector/reference/AutoXAI.md)
  : Auto XAI: Automatic Explainability for Biomarker Models
- [`BayesianTuner`](https://biostat.umed.pl/OmicSelector/reference/BayesianTuner.md)
  : Bayesian Hyperparameter Optimization for Omics
- [`BenchmarkService`](https://biostat.umed.pl/OmicSelector/reference/BenchmarkService.md)
  : BenchmarkService: Nested Cross-Validation with Zero Leakage
- [`DeepLearners`](https://biostat.umed.pl/OmicSelector/reference/DeepLearners.md)
  : Deep Learning Learners for Omics Data
- [`FilterGOF`](https://biostat.umed.pl/OmicSelector/reference/FilterGOF.md)
  : Goodness-of-Fit Filters for Sparse Omics Data
- [`FilterGOF_KS`](https://biostat.umed.pl/OmicSelector/reference/FilterGOF_KS.md)
  : Kolmogorov-Smirnov GOF Filter
- [`FilterHurdle`](https://biostat.umed.pl/OmicSelector/reference/FilterHurdle.md)
  : Hurdle Filter for Zero-Inflated Data
- [`FilterZeroProp`](https://biostat.umed.pl/OmicSelector/reference/FilterZeroProp.md)
  : Zero-Proportion Filter
- [`FrozenComBat`](https://biostat.umed.pl/OmicSelector/reference/FrozenComBat.md)
  : FrozenComBat R6 Class
- [`OmicModalitySpec`](https://biostat.umed.pl/OmicSelector/reference/OmicModalitySpec.md)
  : OmicModalitySpec R6 Class
- [`OmicPipeline`](https://biostat.umed.pl/OmicSelector/reference/OmicPipeline.md)
  : OmicPipeline: Zero-Leakage Feature Selection Pipeline
- [`OmicSelector-package`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector-package.md)
  [`OmicSelector`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector-package.md)
  : OmicSelector: Zero-Leakage Biomarker Discovery Toolkit
- [`OmicSelector_OmicSelector_wrapper()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_OmicSelector_wrapper.md)
  : OmicSelector_OmicSelector (DEPRECATED)
- [`OmicSelector_PCA()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_PCA.md)
  : OmicSelector_PCA
- [`OmicSelector_PCA_3D()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_PCA_3D.md)
  : OmicSelector_PCA_3D
- [`OmicSelector_benchmark()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_benchmark.md)
  : OmicSelector_benchmark
- [`OmicSelector_best_signature_de()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_best_signature_de.md)
  : OmicSelector_best_signature_de
- [`OmicSelector_best_signature_proposals()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_best_signature_proposals.md)
  : OmicSelector_best_signature_proposals
- [`OmicSelector_best_signature_proposals_meta11()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_best_signature_proposals_meta11.md)
  : OmicSelector_best_signature_proposals_meta11
- [`OmicSelector_combat()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_combat.md)
  : OmicSelector_combat
- [`OmicSelector_correct_miRNA_names()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_correct_miRNA_names.md)
  : OmicSelector_correct_miRNA_names
- [`OmicSelector_correlation_plot()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_correlation_plot.md)
  : OmicSelector_correlation_plot
- [`OmicSelector_counts_to_log10tpm()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_counts_to_log10tpm.md)
  : OmicSelector_counts_to_log10tpm
- [`OmicSelector_create_formula()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_create_formula.md)
  : OmicSelector_create_formula
- [`OmicSelector_differential_expression_ttest()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_differential_expression_ttest.md)
  : OmicSelector_differential_expression_ttest
- [`OmicSelector_diverge_color()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_diverge_color.md)
  : OmicSelector_diverge_color
- [`OmicSelector_docker.update_progress()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_docker.update_progress.md)
  : OmicSelector_docker.update_progress
- [`OmicSelector_download_tissue_miRNA_data_from_TCGA()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_download_tissue_miRNA_data_from_TCGA.md)
  : OmicSelector_download_tissue_miRNA_data_from_TCGA
- [`OmicSelector_get_benchmark()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_get_benchmark.md)
  : OmicSelector_get_benchmark
- [`OmicSelector_get_benchmark_methods()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_get_benchmark_methods.md)
  : OmicSelector_get_benchmark_methods
- [`OmicSelector_get_features_from_benchmark()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_get_features_from_benchmark.md)
  : OmicSelector_get_features_from_benchmark
- [`OmicSelector_heatmap.3()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_heatmap.3.md)
  : OmicSelector_heatmap.3
- [`OmicSelector_heatmap()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_heatmap.md)
  : OmicSelector_heatmap
- [`OmicSelector_iteratedRFE()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_iteratedRFE.md)
  : OmicSelector_iteratedRFE
- [`OmicSelector_iteratedRFE_deprecated()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_iteratedRFE_deprecated.md)
  : Wrapper for deprecated OmicSelector_iteratedRFE
- [`OmicSelector_load_datamix()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_load_datamix.md)
  : OmicSelector_load_datamix
- [`OmicSelector_load_extension()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_load_extension.md)
  : OmicSelector_load_extension
- [`OmicSelector_log()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_log.md)
  : OmicSelector_log
- [`OmicSelector_merge_formulas()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_merge_formulas.md)
  : OmicSelector_merge_formulas
- [`OmicSelector_myclust()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_myclust.md)
  : OmicSelector_myclust
- [`OmicSelector_mydist()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_mydist.md)
  : OmicSelector_mydist
- [`OmicSelector_prepare_split()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_prepare_split.md)
  : OmicSelector_prepare_split
- [`OmicSelector_process_tissue_miRNA_TCGA()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_process_tissue_miRNA_TCGA.md)
  : OmicSelector_process_tissue_miRNA_TCGA
- [`OmicSelector_profileplot()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_profileplot.md)
  : OmicSelector_profileplot
- [`OmicSelector_propensity_score_matching()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_propensity_score_matching.md)
  : OmicSelector_propensity_score_matching
- [`OmicSelector_setup()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_setup.md)
  : OmicSelector_setup
- [`OmicSelector_signature_overlap()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_signature_overlap.md)
  : OmicSelector_signature_overlap
- [`OmicSelector_table()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_table.md)
  : OmicSelector_table
- [`OmicSelector_vulcano_plot()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_vulcano_plot.md)
  : OmicSelector_vulcano_plot
- [`OmicSelector_xgboost()`](https://biostat.umed.pl/OmicSelector/reference/OmicSelector_xgboost.md)
  : OmicSelector_xgboost
- [`OmicStackedEnsemble`](https://biostat.umed.pl/OmicSelector/reference/OmicStackedEnsemble.md)
  : OmicStackedEnsemble R6 Class
- [`OmicWeightedEnsemble`](https://biostat.umed.pl/OmicSelector/reference/OmicWeightedEnsemble.md)
  : Simple Weighted Averaging Ensemble
- [`SequentialSelector`](https://biostat.umed.pl/OmicSelector/reference/SequentialSelector.md)
  : Hybrid Sequential Feature Selection (HSFS)
- [`StabilityEnsemble`](https://biostat.umed.pl/OmicSelector/reference/StabilityEnsemble.md)
  : Stability-Based Ensemble Selection
- [`SyntheticData`](https://biostat.umed.pl/OmicSelector/reference/SyntheticData.md)
  : Synthetic Data Generation for Omics
- [`apply_frozen_combat_cv()`](https://biostat.umed.pl/OmicSelector/reference/apply_frozen_combat_cv.md)
  : Apply Frozen ComBat Within Cross-Validation Folds
- [`balance_classes()`](https://biostat.umed.pl/OmicSelector/reference/balance_classes.md)
  : Create Balanced Training Set
- [`build_correlation_adjacency()`](https://biostat.umed.pl/OmicSelector/reference/build_correlation_adjacency.md)
  : Create Correlation-Based Adjacency for GNN
- [`cache`](https://biostat.umed.pl/OmicSelector/reference/cache.md) :
  Split-Aware Caching for OmicSelector 2.0
- [`cache_stats()`](https://biostat.umed.pl/OmicSelector/reference/cache_stats.md)
  : Get Cache Statistics
- [`cached_filter()`](https://biostat.umed.pl/OmicSelector/reference/cached_filter.md)
  : Cached Filter Computation
- [`calibration`](https://biostat.umed.pl/OmicSelector/reference/calibration.md)
  : Calibration Metrics for OmicSelector 2.0
- [`calibration_summary()`](https://biostat.umed.pl/OmicSelector/reference/calibration_summary.md)
  : Calibration Summary for Model Results
- [`check_batch_correction_leakage()`](https://biostat.umed.pl/OmicSelector/reference/check_batch_correction_leakage.md)
  : Check for Batch Correction Leakage
- [`check_dl_availability()`](https://biostat.umed.pl/OmicSelector/reference/check_dl_availability.md)
  : Check Deep Learning Availability
- [`check_feature_correlations()`](https://biostat.umed.pl/OmicSelector/reference/check_feature_correlations.md)
  : Check Feature Correlations
- [`clear_cache()`](https://biostat.umed.pl/OmicSelector/reference/clear_cache.md)
  : Clear OmicSelector Cache
- [`compare_gof_filters()`](https://biostat.umed.pl/OmicSelector/reference/compare_gof_filters.md)
  : Compare GOF Filters on Task
- [`compute_ece()`](https://biostat.umed.pl/OmicSelector/reference/compute_ece.md)
  : Compute Expected Calibration Error (ECE)
- [`compute_nogueira_stability()`](https://biostat.umed.pl/OmicSelector/reference/compute_nogueira_stability.md)
  : Compute Nogueira Stability Index
- [`compute_shap_with_warnings()`](https://biostat.umed.pl/OmicSelector/reference/compute_shap_with_warnings.md)
  : Compute SHAP Values with Correlation Warnings
- [`compute_stability_from_resample()`](https://biostat.umed.pl/OmicSelector/reference/compute_stability_from_resample.md)
  : Compute Stability from ResampleResult
- [`create_explainer()`](https://biostat.umed.pl/OmicSelector/reference/create_explainer.md)
  : Create Model Explainer
- [`create_frozen_combat_pipeop()`](https://biostat.umed.pl/OmicSelector/reference/create_frozen_combat_pipeop.md)
  : Create a Frozen ComBat PipeOp for mlr3pipelines
- [`create_hsfs_selector()`](https://biostat.umed.pl/OmicSelector/reference/create_hsfs_selector.md)
  : Create a Hybrid Sequential Feature Selector
- [`create_mlp_learner()`](https://biostat.umed.pl/OmicSelector/reference/create_mlp_learner.md)
  : Create MLP Learner via mlr3torch
- [`create_omic_cache()`](https://biostat.umed.pl/OmicSelector/reference/create_omic_cache.md)
  : Create Split-Aware Cache
- [`create_report_data()`](https://biostat.umed.pl/OmicSelector/reference/create_report_data.md)
  : Create Report Data Schema
- [`create_stability_ensemble()`](https://biostat.umed.pl/OmicSelector/reference/create_stability_ensemble.md)
  : Create a Stability Ensemble with Presets
- [`decompose_brier()`](https://biostat.umed.pl/OmicSelector/reference/decompose_brier.md)
  : Decompose Brier Score
- [`deprecated`](https://biostat.umed.pl/OmicSelector/reference/deprecated.md)
  : Deprecated and Leaky Functions
- [`export_bundle()`](https://biostat.umed.pl/OmicSelector/reference/export_bundle.md)
  : Create Complete Export Bundle
- [`export_onnx()`](https://biostat.umed.pl/OmicSelector/reference/export_onnx.md)
  : Export Model to ONNX Format
- [`export_vetiver()`](https://biostat.umed.pl/OmicSelector/reference/export_vetiver.md)
  : Export Model as Vetiver for Deployment
- [`extract_features_from_resample()`](https://biostat.umed.pl/OmicSelector/reference/extract_features_from_resample.md)
  : Extract Features from All Folds in ResampleResult
- [`extract_selected_features()`](https://biostat.umed.pl/OmicSelector/reference/extract_selected_features.md)
  : Extract Selected Features from Trained GraphLearner
- [`feature_importance()`](https://biostat.umed.pl/OmicSelector/reference/feature_importance.md)
  : Compute Permutation Feature Importance
- [`fit_isotonic_calibration()`](https://biostat.umed.pl/OmicSelector/reference/fit_isotonic_calibration.md)
  : Isotonic Regression Calibration
- [`fit_platt_scaling()`](https://biostat.umed.pl/OmicSelector/reference/fit_platt_scaling.md)
  : Platt Scaling (Logistic Calibration)
- [`fit_temperature_scaling()`](https://biostat.umed.pl/OmicSelector/reference/fit_temperature_scaling.md)
  : Temperature Scaling
- [`frozen-combat`](https://biostat.umed.pl/OmicSelector/reference/frozen-combat.md)
  : Frozen ComBat for Leakage-Free Batch Correction
- [`frozen_combat_correct()`](https://biostat.umed.pl/OmicSelector/reference/frozen_combat_correct.md)
  : Convenience function for frozen ComBat correction
- [`generate_cache_key()`](https://biostat.umed.pl/OmicSelector/reference/generate_cache_key.md)
  : Generate Split-Aware Cache Key
- [`generate_tripod_report()`](https://biostat.umed.pl/OmicSelector/reference/generate_tripod_report.md)
  : Generate TRIPOD+AI Report
- [`get_consensus_features()`](https://biostat.umed.pl/OmicSelector/reference/get_consensus_features.md)
  : Get Consensus Features from Best Signature
- [`get_modality_info()`](https://biostat.umed.pl/OmicSelector/reference/get_modality_info.md)
  : Get Modality Information
- [`get_optimal_params()`](https://biostat.umed.pl/OmicSelector/reference/get_optimal_params.md)
  : Get Optimal Hyperparameters from AutoTuner
- [`get_parallel_status()`](https://biostat.umed.pl/OmicSelector/reference/get_parallel_status.md)
  : Get Current Parallelization Status
- [`get_reliable_shap_features()`](https://biostat.umed.pl/OmicSelector/reference/get_reliable_shap_features.md)
  : Get Reliable SHAP Features
- [`get_selected_features_per_fold()`](https://biostat.umed.pl/OmicSelector/reference/get_selected_features_per_fold.md)
  : Get Selected Features Per Fold
- [`interpretability`](https://biostat.umed.pl/OmicSelector/reference/interpretability.md)
  : Model Interpretability for OmicSelector 2.0
- [`list_deprecated_functions()`](https://biostat.umed.pl/OmicSelector/reference/list_deprecated_functions.md)
  : List all deprecated functions
- [`load_bundle()`](https://biostat.umed.pl/OmicSelector/reference/load_bundle.md)
  : Load Exported Model Bundle
- [`make_autotuner_glmnet()`](https://biostat.umed.pl/OmicSelector/reference/make_autotuner_glmnet.md)
  : Create Bayesian-Optimized AutoTuner for glmnet
- [`make_autotuner_lightgbm()`](https://biostat.umed.pl/OmicSelector/reference/make_autotuner_lightgbm.md)
  : Create Bayesian-Optimized AutoTuner for LightGBM
- [`make_autotuner_ranger()`](https://biostat.umed.pl/OmicSelector/reference/make_autotuner_ranger.md)
  : Create Bayesian-Optimized AutoTuner for Random Forest
- [`make_autotuner_xgboost()`](https://biostat.umed.pl/OmicSelector/reference/make_autotuner_xgboost.md)
  : Create Bayesian-Optimized AutoTuner for XGBoost
- [`make_gnn_learner()`](https://biostat.umed.pl/OmicSelector/reference/make_gnn_learner.md)
  : Create GNN Learner for Pathway-Aware Classification
- [`make_gof_filter()`](https://biostat.umed.pl/OmicSelector/reference/make_gof_filter.md)
  : Create GOF Filter
- [`make_mlp_learner()`](https://biostat.umed.pl/OmicSelector/reference/make_mlp_learner.md)
  : Create Omics-Optimized MLP Learner
- [`make_tabtransformer_learner()`](https://biostat.umed.pl/OmicSelector/reference/make_tabtransformer_learner.md)
  : Create TabTransformer Learner
- [`memoize_with_split()`](https://biostat.umed.pl/OmicSelector/reference/memoize_with_split.md)
  : Memoize Function with Split Context
- [`merge_omics_data()`](https://biostat.umed.pl/OmicSelector/reference/merge_omics_data.md)
  : Merge Multi-Omics Data for Analysis
- [`model-export`](https://biostat.umed.pl/OmicSelector/reference/model-export.md)
  : Model Export for OmicSelector 2.0
- [`multi-omics`](https://biostat.umed.pl/OmicSelector/reference/multi-omics.md)
  : Multi-Omics Support for OmicSelector 2.0
- [`noise_augment()`](https://biostat.umed.pl/OmicSelector/reference/noise_augment.md)
  : Gaussian Noise Augmentation
- [`omic_benchmark()`](https://biostat.umed.pl/OmicSelector/reference/omic_benchmark.md)
  : Create benchmark service from OmicPipeline
- [`omic_pipeline()`](https://biostat.umed.pl/OmicSelector/reference/omic_pipeline.md)
  : Quick pipeline creation from data
- [`parallel`](https://biostat.umed.pl/OmicSelector/reference/parallel.md)
  : Parallelization Support for OmicSelector 2.0
- [`partial_dependence()`](https://biostat.umed.pl/OmicSelector/reference/partial_dependence.md)
  : Compute Partial Dependence
- [`plot_signature_tradeoffs()`](https://biostat.umed.pl/OmicSelector/reference/plot_signature_tradeoffs.md)
  : Plot Signature Selection Trade-offs
- [`plot_xai_importance()`](https://biostat.umed.pl/OmicSelector/reference/plot_xai_importance.md)
  : Plot XAI Feature Importance
- [`print(`*`<CalibrationResult>`*`)`](https://biostat.umed.pl/OmicSelector/reference/print.CalibrationResult.md)
  : Print method for CalibrationResult
- [`print(`*`<CorrelationCheck>`*`)`](https://biostat.umed.pl/OmicSelector/reference/print.CorrelationCheck.md)
  : Print Correlation Check
- [`print(`*`<FeatureImportance>`*`)`](https://biostat.umed.pl/OmicSelector/reference/print.FeatureImportance.md)
  : Print Feature Importance
- [`print(`*`<LeakageValidation>`*`)`](https://biostat.umed.pl/OmicSelector/reference/print.LeakageValidation.md)
  : Print method for LeakageValidation
- [`print(`*`<NestedCVResult>`*`)`](https://biostat.umed.pl/OmicSelector/reference/print.NestedCVResult.md)
  : Print method for NestedCVResult
- [`print(`*`<NogueiraStability>`*`)`](https://biostat.umed.pl/OmicSelector/reference/print.NogueiraStability.md)
  : Print method for NogueiraStability
- [`print(`*`<OmicBenchmarkResult>`*`)`](https://biostat.umed.pl/OmicSelector/reference/print.OmicBenchmarkResult.md)
  : Print method for OmicBenchmarkResult
- [`print(`*`<OmicsInput>`*`)`](https://biostat.umed.pl/OmicSelector/reference/print.OmicsInput.md)
  : Print method for OmicsInput
- [`print(`*`<ReportData>`*`)`](https://biostat.umed.pl/OmicSelector/reference/print.ReportData.md)
  : Print method for ReportData
- [`print(`*`<SignatureSelectionResult>`*`)`](https://biostat.umed.pl/OmicSelector/reference/print.SignatureSelectionResult.md)
  : Print Method for Signature Selection Result
- [`print_shap_warnings()`](https://biostat.umed.pl/OmicSelector/reference/print_shap_warnings.md)
  : Print SHAP Warnings Report
- [`print_xai_summary()`](https://biostat.umed.pl/OmicSelector/reference/print_xai_summary.md)
  : Print XAI Summary
- [`register_gof_filters()`](https://biostat.umed.pl/OmicSelector/reference/register_gof_filters.md)
  : Register GOF Filters in mlr3
- [`reliability_diagram_data()`](https://biostat.umed.pl/OmicSelector/reference/reliability_diagram_data.md)
  : Create Reliability Diagram Data
- [`report`](https://biostat.umed.pl/OmicSelector/reference/report.md) :
  TRIPOD+AI Report Generation for OmicSelector 2.0
- [`reset_parallel()`](https://biostat.umed.pl/OmicSelector/reference/reset_parallel.md)
  : Reset Parallelization to Sequential
- [`run_bayesian_benchmark()`](https://biostat.umed.pl/OmicSelector/reference/run_bayesian_benchmark.md)
  : Run Bayesian Optimization Benchmark
- [`run_dl_benchmark()`](https://biostat.umed.pl/OmicSelector/reference/run_dl_benchmark.md)
  : Deep Learning Benchmark
- [`select_best_signature()`](https://biostat.umed.pl/OmicSelector/reference/select_best_signature.md)
  : Select Best Biomarker Signature from Nested CV Results
- [`setup_parallel()`](https://biostat.umed.pl/OmicSelector/reference/setup_parallel.md)
  : Configure Parallelization for OmicSelector
- [`shap_values()`](https://biostat.umed.pl/OmicSelector/reference/shap_values.md)
  : Compute SHAP-like Values
- [`shap_warnings`](https://biostat.umed.pl/OmicSelector/reference/shap_warnings.md)
  : Correlation-Aware SHAP Interpretation
- [`signature-selection`](https://biostat.umed.pl/OmicSelector/reference/signature-selection.md)
  : Signature Selection: Multi-Objective Best Biomarker Selection
- [`smote_augment()`](https://biostat.umed.pl/OmicSelector/reference/smote_augment.md)
  : SMOTE Augmentation for Omics Data
- [`stability`](https://biostat.umed.pl/OmicSelector/reference/stability.md)
  : Nogueira Stability Index for Feature Selection
- [`stack_omics()`](https://biostat.umed.pl/OmicSelector/reference/stack_omics.md)
  : Create Multi-Omics Stacked Ensemble (Convenience Function)
- [`tabddpm_generate()`](https://biostat.umed.pl/OmicSelector/reference/tabddpm_generate.md)
  : TabDDPM Synthetic Data Generator
- [`torch-learners`](https://biostat.umed.pl/OmicSelector/reference/torch-learners.md)
  : mlr3torch Learner Integration for OmicSelector 2.0
- [`validate_no_leakage()`](https://biostat.umed.pl/OmicSelector/reference/validate_no_leakage.md)
  : Validate pipeline for leakage risks
- [`validate_omics_input()`](https://biostat.umed.pl/OmicSelector/reference/validate_omics_input.md)
  : Validate Multi-Omics Input
- [`validate_synthetic()`](https://biostat.umed.pl/OmicSelector/reference/validate_synthetic.md)
  : Validate Synthetic Data Quality
- [`with_parallel()`](https://biostat.umed.pl/OmicSelector/reference/with_parallel.md)
  : With Parallel Scope
- [`xai_correlations()`](https://biostat.umed.pl/OmicSelector/reference/xai_correlations.md)
  : Compute Correlation Diagnostics for Features
- [`xai_explainer_mlr3()`](https://biostat.umed.pl/OmicSelector/reference/xai_explainer_mlr3.md)
  : Create DALEX Explainer from mlr3 Learner
- [`xai_importance()`](https://biostat.umed.pl/OmicSelector/reference/xai_importance.md)
  : Compute Permutation Feature Importance
- [`xai_pdp()`](https://biostat.umed.pl/OmicSelector/reference/xai_pdp.md)
  : Compute Partial Dependence Plots
- [`xai_pipeline()`](https://biostat.umed.pl/OmicSelector/reference/xai_pipeline.md)
  : Run Complete XAI Pipeline
- [`xai_shap()`](https://biostat.umed.pl/OmicSelector/reference/xai_shap.md)
  : Compute SHAP Values for Observations
