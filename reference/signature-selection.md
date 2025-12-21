# Signature Selection: Multi-Objective Best Biomarker Selection

Functions for selecting the best biomarker signature from nested CV
results. Implements multiple selection strategies: constrained 1SE rule,
weighted scoring, and Pareto optimality for transparent multi-objective
trade-offs.

## Details

\## Clinical Context

In biomarker discovery, selecting the "best" signature requires
balancing: - \*\*Performance\*\* (AUC, accuracy): Predictive power -
\*\*Stability\*\*: Reproducibility across resampling (Nogueira Index) -
\*\*Parsimony\*\*: Fewer features = cheaper, more practical clinical
panels

This module provides three selection modes: 1. \`constrained_1se\`:
Conservative selection within 1 standard error of best 2. \`weighted\`:
Explicit weighted combination of objectives 3. \`pareto\`: Return
non-dominated solutions for human inspection
