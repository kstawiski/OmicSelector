# Hemolysis-Aware Corrections for Biomarker Panels

Utilities for simple deployment-facing hemolysis mitigations used in
biomarker-panel robustness benchmarks:

- a Blondal-style proxy pre-filter tuned on training controls only;

- a weighted centered log-ratio (CLR) transform that down-weights
  hemolysis-sensitive features while preserving within-sample centering.

The default proxy uses the eLife-14 panel overlap markers \`miR-92a-3p /
miR-320d\` because canonical \`miR-23a / miR-451a\` are not available in
the frozen 14-miRNA panel.
