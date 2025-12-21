# Frozen ComBat for Leakage-Free Batch Correction

Implements a "frozen" version of ComBat batch correction that prevents
data leakage during external validation. Standard ComBat estimates batch
parameters from ALL data, which causes leakage when applied to
train+test together. Frozen ComBat separates parameter estimation
(training only) from application (any data).

## Details

\## Why Frozen ComBat?

Standard ComBat workflow: “\` combined_data \<- rbind(train, test)
corrected \<- ComBat(combined_data, batch) \# LEAKAGE: test data
influences correction “\`

Frozen ComBat workflow: “\` frozen \<- FrozenComBat\$new()
frozen\$fit(train_data, train_batch) \# Learn parameters from training
only corrected_train \<- frozen\$transform(train_data, train_batch)
corrected_test \<- frozen\$transform(test_data, test_batch) \# Apply
frozen parameters “\`

## References

Johnson, W.E., Li, C., Rabinovic, A. (2007). Adjusting batch effects in
microarray expression data using empirical Bayes methods. Biostatistics,
8(1), 118-127.
