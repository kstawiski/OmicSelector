# mlr3torch Learner Integration for OmicSelector 2.0

Factory functions for creating mlr3torch-based deep learning models.
Phase 2 scope: Multi-Layer Perceptron (MLP) only.

## Details

Why MLP only for Phase 2: - Tabular omics data typically doesn't have
spatial/temporal ordering - CNNs assume local structure (adjacent
features correlated) - misleading for most omics - MLPs are the
appropriate baseline for unordered feature vectors - Focus is on proper
nested CV and leakage prevention, not architecture variety
