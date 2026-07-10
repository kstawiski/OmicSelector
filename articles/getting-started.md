# Getting Started with OmicSelector

## Overview

This vignette shows a full workflow on a **tiny TCGA subset**:

1.  Load data
2.  Create an `OmicPipeline`
3.  Build a leakage-safe `GraphLearner`
4.  Run nested CV
5.  Fit a final model and extract selected features

## Load Package and Data

``` r

library(OmicSelector)

# TCGA miRNA dataset shipped with the package
# (data frame with clinical + miRNA features)
data("original_TCGA_data", package = "OmicSelector")
```

## Prepare a Tiny TCGA Subset

``` r

feature_cols <- grep("^hsa\\.", names(original_TCGA_data), value = TRUE)
feature_cols <- head(feature_cols, 200)

# Keep only patient id, target, and a small feature subset
subset_df <- original_TCGA_data[, c("patient", "sample_type", feature_cols), drop = FALSE]
subset_df <- as.data.frame(subset_df)
subset_df$sample_type <- factor(subset_df$sample_type)

# Balance classes (small and fast)
set.seed(1)
idx_tumor <- which(subset_df$sample_type == "PrimaryTumor")
idx_normal <- which(subset_df$sample_type == "SolidTissueNormal")
subset_df <- subset_df[c(sample(idx_tumor, 50), sample(idx_normal, 50)), ]
subset_df <- subset_df[sample(nrow(subset_df)), ]
```

## Create Pipeline

``` r

pipeline <- OmicPipeline$new(
  data = subset_df,
  target = "sample_type",
  positive = "PrimaryTumor",
  patient_id = "patient"
)
```

## Build GraphLearner (with Screening)

``` r

learner <- pipeline$create_graph_learner(
  filter = "anova",
  model = "rpart",
  n_features = 20,
  screening = TRUE,
  screening_frac = 0.2
)
```

## Optional: Deep Learning + Autoencoder

``` r

# Requires torch + mlr3torch
dl_learner <- pipeline$create_graph_learner(
  filter = "anova",
  model = "mlp",
  n_features = 50,
  autoencoder = list(
    latent_dim = 32,
    hidden_layers = c(128, 64),
    epochs = 50,
    batch_size = 64,
    early_stopping = TRUE,
    patience = 10
  )
)
```

## Run Nested CV

``` r

result <- pipeline$benchmark(
  learners = learner,
  outer_folds = 3,
  inner_folds = 2,
  seed = 42,
  parallel = TRUE,
  threads = 1,
  cache_dir = "cache"
)

print(result)
```

## Fit Final Model and Extract Features

``` r

fit <- pipeline$fit(
  learner = learner,
  seed = 42,
  threads = 1
)

fit$selected_features
```

## Next Steps

- Increase folds for a more reliable estimate
- Swap filters (`"mrmr"`, `"relief"`, `"correlation"`)
- Use multi‑omics input by passing a **named list** of data.frames
