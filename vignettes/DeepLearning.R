## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, message=FALSE, warning=FALSE,
  comment = "#>"
)
knitr::opts_chunk$set(fig.width=12, fig.height=8)
knitr::opts_chunk$set(tidy.opts=list(width.cutoff=150),tidy=TRUE)
options(rgl.useNULL = TRUE)
options(warn=-1)
suppressMessages(library(dplyr))
set.seed(1)
options(knitr.table.format = "html")
library(OmicSelector)

## -----------------------------------------------------------------------------
OmicSelector::OmicSelector_load_extension("deeplearning")

## ----eval=F-------------------------------------------------------------------
#  OmicSelector_deep_learning = function(selected_miRNAs = ".",
#                                        wd = getwd(),
#                                        SMOTE = F,
#                                        keras_batch_size = 64,
#                                        clean_temp_files = T,
#                                        save_threshold_trainacc = 0.85,
#                                        save_threshold_testacc = 0.8,
#                                        keras_epochae = 5000,
#                                        keras_epoch = 2000,
#                                        keras_patience = 50,
#                                        hyperparameters = expand.grid(
#                                          layer1 = seq(3, 11, by = 2),
#                                          layer2 = c(0, seq(3, 11, by = 2)),
#                                          layer3 c(0, seq(3, 11, by = 2)),
#                                          activation_function_layer1 = c("relu", "sigmoid"),
#                                          activation_function_layer2 = c("relu", "sigmoid"),
#                                          activation_function_layer3 = c("relu", "sigmoid"),
#                                          dropout_layer1 = c(0, 0.1),
#                                          dropout_layer2 = c(0, 0.1),
#                                          dropout_layer3 = c(0),
#                                          layer1_regularizer = c(T, F),
#                                          layer2_regularizer = c(T, F),
#                                          layer3_regularizer = c(T, F),
#                                          optimizer = c("adam", "rmsprop", "sgd"),
#                                          autoencoder = c(0, 7, -7),
#                                          balanced = SMOTE,
#                                          formula = as.character(OmicSelector_create_formula(selected_miRNAs))[3],
#                                          scaled = c(T, F),
#                                          stringsAsFactors = F
#                                        ),
#                                        add_features_to_predictions = F,
#                                        keras_threads = ceiling(parallel::detectCores() /
#                                                                  2),
#                                        start = 1,
#                                        end = nrow(hyperparameters),
#                                        output_file = "deeplearning_results.csv",
#                                        save_all_vars = F)
#  OmicSelector_load_datamix()

## ----eval=T-------------------------------------------------------------------
balanced = F # not balanced
selected_miRNAs = c("hsa.a","hsa.b","hsa.c") # selected features
hyperparameters = expand.grid(
  layer1 = seq(2, 10, by = 1),
  layer2 = c(0),
  layer3 = c(0),
  activation_function_layer1 = c("relu", "sigmoid", "selu"),
  activation_function_layer2 = c("relu"),
  activation_function_layer3 = c("relu"),
  dropout_layer1 = c(0, 0.1),
  dropout_layer2 = c(0),
  dropout_layer3 = c(0),
  layer1_regularizer = c(T, F),
  layer2_regularizer = c(F),
  layer3_regularizer = c(F),
  optimizer = c("adam", "rmsprop", "sgd"),
  autoencoder = c(0, -7, 7),
  balanced = balanced,
  formula = as.character(OmicSelector_create_formula(selected_miRNAs))[3],
  scaled = c(T, F),
  stringsAsFactors = F
)
# DT::datatable(hyperparameters, 
#          extensions = c('FixedColumns',"FixedHeader"),
#           options = list(scrollX = TRUE, 
#                          paging=FALSE,
#                          fixedHeader=TRUE))

## ----eval=F-------------------------------------------------------------------
#  hyperparameters_part1 = expand.grid(
#    layer1 = seq(2, 10, by = 1),
#    layer2 = c(0),
#    layer3 = c(0),
#    activation_function_layer1 = c("relu", "sigmoid", "selu"),
#    activation_function_layer2 = c("relu"),
#    activation_function_layer3 = c("relu"),
#    dropout_layer1 = c(0, 0.1),
#    dropout_layer2 = c(0),
#    dropout_layer3 = c(0),
#    layer1_regularizer = c(T, F),
#    layer2_regularizer = c(F),
#    layer3_regularizer = c(F),
#    optimizer = c("adam", "rmsprop", "sgd"),
#    autoencoder = c(0),
#    balanced = balanced,
#    formula = as.character(OmicSelector_create_formula(selected_miRNAs))[3],
#    scaled = c(T, F),
#    stringsAsFactors = F
#  )
#  hyperparameters_part2 = expand.grid(
#    layer1 = seq(3, 11, by = 2),
#    layer2 = c(seq(3, 11, by = 2)),
#    layer3 = c(seq(0, 11, by = 2)),
#    activation_function_layer1 = c("relu", "sigmoid", "selu"),
#    activation_function_layer2 = c("relu", "sigmoid", "selu"),
#    activation_function_layer3 = c("relu", "sigmoid", "selu"),
#    dropout_layer1 = c(0, 0.1),
#    dropout_layer2 = c(0),
#    dropout_layer3 = c(0),
#    layer1_regularizer = c(T, F),
#    layer2_regularizer = c(F),
#    layer3_regularizer = c(F),
#    optimizer = c("adam", "rmsprop", "sgd"),
#    autoencoder = c(0),
#    balanced = balanced,
#    formula = as.character(OmicSelector_create_formula(selected_miRNAs))[3],
#    scaled = c(T, F),
#    stringsAsFactors = F
#  )
#  hyperparameters = rbind(hyperparameters_part1, hyperparameters_part2)

## ----eval=F-------------------------------------------------------------------
#  hyperparameters_part1 = expand.grid(
#    layer1 = seq(2, 10, by = 1),
#    layer2 = c(0),
#    layer3 = c(0),
#    activation_function_layer1 = c("relu", "sigmoid", "selu"),
#    activation_function_layer2 = c("relu"),
#    activation_function_layer3 = c("relu"),
#    dropout_layer1 = c(0, 0.1),
#    dropout_layer2 = c(0),
#    dropout_layer3 = c(0),
#    layer1_regularizer = c(T, F),
#    layer2_regularizer = c(F),
#    layer3_regularizer = c(F),
#    optimizer = c("adam", "rmsprop", "sgd"),
#    autoencoder = c(0, -7, 7),
#    balanced = balanced,
#    formula = as.character(OmicSelector_create_formula(selected_miRNAs))[3],
#    scaled = c(T, F),
#    stringsAsFactors = F
#  )
#  hyperparameters_part2 = expand.grid(
#    layer1 = seq(3, 11, by = 2),
#    layer2 = c(seq(3, 11, by = 2)),
#    layer3 = c(seq(0, 11, by = 2)),
#    activation_function_layer1 = c("relu", "sigmoid", "selu"),
#    activation_function_layer2 = c("relu", "sigmoid", "selu"),
#    activation_function_layer3 = c("relu", "sigmoid", "selu"),
#    dropout_layer1 = c(0, 0.1),
#    dropout_layer2 = c(0),
#    dropout_layer3 = c(0),
#    layer1_regularizer = c(T, F),
#    layer2_regularizer = c(F),
#    layer3_regularizer = c(F),
#    optimizer = c("adam", "rmsprop", "sgd"),
#    autoencoder = c(0, -7, 7),
#    balanced = balanced,
#    formula = as.character(OmicSelector_create_formula(selected_miRNAs))[3],
#    scaled = c(T, F),
#    stringsAsFactors = F
#  )
#  hyperparameters = rbind(hyperparameters_part1, hyperparameters_part2)

## -----------------------------------------------------------------------------
data("orginal_TCGA_data")
suppressWarnings(suppressMessages(library(dplyr)))
cancer_cases = filter(orginal_TCGA_data, primary_site == "Pancreas" & sample_type == "PrimaryTumor")
control_cases = filter(orginal_TCGA_data, sample_type == "SolidTissueNormal")
cancer_cases$Class = "Case"
control_cases$Class = "Control"
dataset = rbind(cancer_cases, control_cases)
# Not run:
# DE = OmicSelector_differential_expression_ttest(ttpm_features = dplyr::select(dataset, starts_with("hsa")), classes = dataset$Class, mode = "logtpm")
# significant = DE$miR[DE$`p-value Bonferroni`<0.05]
selected_miRNAs = make.names(c('hsa-miR-192-5p', 'hsa-let-7g-5p', 'hsa-let-7a-5p', 'hsa-miR-194-5p', 'hsa-miR-122-5p', 'hsa-miR-340-5p', 'hsa-miR-26b-5p')) # some selected miRNAs
dataset = dataset[sample(1:nrow(dataset),200),] # sample 100 random cases to make it quicker
OmicSelector_table(table(dataset$Class), col.names = c("Class", "Number of cases"))
merged = OmicSelector_prepare_split(metadane = dplyr::select(dataset, -starts_with("hsa")), ttpm = dplyr::select(dataset, selected_miRNAs))
knitr::kable(table(merged$Class, merged$mix))

## ----include = FALSE----------------------------------------------------------
OmicSelector_load_datamix()

## ----eval=FALSE---------------------------------------------------------------
#  SMOTE = F # use unbalanced set
#  deep_learning_results = OmicSelector_deep_learning(selected_miRNAs = selected_miRNAs, start = 5, end = 10) # use default set of hyperparameters and options
#  DT::datatable(deep_learning_results,
#            extensions = c('FixedColumns',"FixedHeader"),
#             options = list(scrollX = TRUE,
#                            paging=FALSE,
#                            fixedHeader=TRUE))

## ----echo = FALSE-------------------------------------------------------------
deep_learning_results = data.table::fread("deeplearning_results.cs")
DT::datatable(deep_learning_results, 
          extensions = c('FixedColumns',"FixedHeader"),
           options = list(scrollX = TRUE, 
                          paging=FALSE,
                          fixedHeader=TRUE))

## ----eval=F-------------------------------------------------------------------
#  OmicSelector_deep_learning_predict = function(model_path = "our_models/model5.zip",
#                                                new_dataset = data.table::fread("Data/ks_data.csv"),
#                                                new_scaling = F,
#                                                old_train_csv_to_restore_scaling = NULL,
#                                                override_cutoff = NULL,
#                                                blinded = F)

