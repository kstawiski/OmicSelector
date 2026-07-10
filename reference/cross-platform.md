# Cross-Platform Transfer Learning Utilities

Utilities for adapting feature distributions between measurement
platforms before model transfer. This is designed for settings where
standard frozen batch correction is not applicable because the source
training set contains only one platform level.

## Details

\`CrossPlatformAdapter\` implements lightweight platform adaptation
strategies that can be applied before training on a source platform and
predicting on a target platform:

\- \`"rank"\`: within-sample percentile ranks across features -
\`"quantile"\`: quantile normalization to the source reference
distribution - \`"zscore"\`: within-sample z-scoring - \`"reference"\`:
subtraction of per-sample reference miRNA summary - \`"combat_pooled"\`:
pooled ComBat correction across source and target

The pooled ComBat strategy is semi-supervised because it needs some
target labels to preserve biology while correcting platform effects.
