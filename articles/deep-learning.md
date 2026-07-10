# Deep Learning with OmicSelector

## Overview

This vignette demonstrates a full deep learning workflow with
OmicSelector, including autoencoder feature compression, mlr3torch
training, checkpoint export/import, and fine‑tuning.

**Requirements:** - `torch` - `mlr3torch`

``` r

install.packages("torch")
torch::install_torch()
install.packages("mlr3torch")
```

## Load Package and Data

``` r

library(OmicSelector)

data("original_TCGA_data", package = "OmicSelector")
```

## Prepare a Small TCGA Subset

``` r

feature_cols <- grep("^hsa\\.", names(original_TCGA_data), value = TRUE)
feature_cols <- head(feature_cols, 300)

subset_df <- original_TCGA_data[, c("patient", "sample_type", feature_cols), drop = FALSE]
subset_df <- as.data.frame(subset_df)
subset_df$sample_type <- factor(subset_df$sample_type)

set.seed(1)
idx_tumor <- which(subset_df$sample_type == "PrimaryTumor")
idx_normal <- which(subset_df$sample_type == "SolidTissueNormal")
subset_df <- subset_df[c(sample(idx_tumor, 60), sample(idx_normal, 60)), ]
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

## Build Deep Learning GraphLearner

This pipeline adds a torch autoencoder *before* feature selection and
trains an MLP via mlr3torch.

``` r

learner <- pipeline$create_graph_learner(
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
  threads = 1
)

print(result)
```

## Fit Final Model

``` r

fit <- pipeline$fit(learner, seed = 42)
fit$selected_features
```

## Export / Import Checkpoint

``` r

export_omicfit_checkpoint(fit, "mlp_checkpoint.pt")

# Load into a trained learner (for inference)
fit <- import_omicfit_checkpoint(fit, "mlp_checkpoint.pt")
```

## Fine‑tune on a New Task

``` r

task <- pipeline$get_task()

finetuned <- finetune_mlr3torch_checkpoint(
  learner = learner,
  task = task,
  path = "mlp_checkpoint.pt",
  epochs = 30
)

# Or stay at the OmicFit level
finetuned_fit <- finetune_omicfit_checkpoint(
  fit = fit,
  task = task,
  path = "mlp_checkpoint.pt",
  epochs = 30
)
```

## Notes

- Autoencoders can also be trained separately via
  [`autoencoder_fit()`](https://kstawiski.github.io/OmicSelector/reference/autoencoder_fit.md)
  and passed into the pipeline using `pretrained` and `freeze_encoder`.
- For transfer learning across datasets, export checkpoints and re‑train
  with
  [`finetune_mlr3torch_checkpoint()`](https://kstawiski.github.io/OmicSelector/reference/finetune_mlr3torch_checkpoint.md).
