# Do wywalenia potem: 
studyMirs = make.names(c('hsa-miR-192-5p', 'hsa-let-7g-5p', 'hsa-let-7a-5p', 'hsa-let-7d-5p', 'hsa-miR-194-5p', 'hsa-miR-98-5p', 'hsa-let-7f-5p', 'hsa-miR-122-5p', 'hsa-miR-340-5p', 'hsa-miR-26b-5p'))
norm3_1 = make.names(c('hsa-miR-17-5p', 'hsa-miR-92a-3p', 'hsa-miR-199a-3p'))
norm2 = make.names(c('hsa-miR-28-3p', 'hsa-miR-92a-3p'))
normQ = make.names('hsa-miR-23a-3p')
all_norm = make.names(c('hsa-miR-17-5p', 'hsa-miR-92a-3p', 'hsa-miR-199a-3p', 'hsa-miR-23a-3p', 'hsa-miR-28-3p'))

# Parameters:
balanced = F
nazwa_konfiguracji = "init.csv"
selected_miRNAs = studyMirs

# Code:
options(warn = -1)
if(file.exists("task.log")) { file.remove("task.log") }
con <- file("task.log")
sink(con, append=TRUE)
sink(con, append=TRUE, type = "message")
library(OmicSelector); set.seed(1);
if(!dir.exists("/OmicSelector/temp")) { dir.create("/OmicSelector/temp") }
OmicSelector_load_extension("deeplearning")
library(data.table)

dane = OmicSelector_load_datamix()
if(balanced == F) {
  t = data.table::fread("mixed_train.csv")
} else {
  t = data.table::fread("mixed_train_balanced.csv")
}


# autoencoder ma softmax na deep feature
hyperparameters_part1 = expand.grid(layer1 = seq(2,10, by = 1), layer2 = c(0), layer3 = c(0),
                                    activation_function_layer1 = c("relu","sigmoid","selu"), activation_function_layer2 = c("relu"), activation_function_layer3 = c("relu"),
                                    dropout_layer1 = c(0, 0.1), dropout_layer2 = c(0), dropout_layer3 = c(0),
                                    layer1_regularizer = c(T,F), layer2_regularizer = c(F), layer3_regularizer = c(F),
                                    optimizer = c("adam","rmsprop","sgd"), autoencoder = c(0), balanced = balanced, formula = as.character(OmicSelector_create_formula(selected_miRNAs))[3], scaled = c(T,F),
                                    stringsAsFactors = F)
hyperparameters_part2 = expand.grid(layer1 = seq(3,11, by = 2), layer2 = c(seq(3,11, by = 2)), layer3 = c(seq(0,11, by = 2)),
                                    activation_function_layer1 = c("relu","sigmoid","selu"), activation_function_layer2 = c("relu","sigmoid","selu"), activation_function_layer3 = c("relu","sigmoid","selu"),
                                    dropout_layer1 = c(0, 0.1), dropout_layer2 = c(0), dropout_layer3 = c(0),
                                    layer1_regularizer = c(T,F), layer2_regularizer = c(F), layer3_regularizer = c(F),
                                    optimizer = c("adam","rmsprop","sgd"), autoencoder = c(0), balanced = balanced, formula = as.character(OmicSelector_create_formula(selected_miRNAs))[3], scaled = c(T,F),
                                    stringsAsFactors = F)
hyperparameters = rbind(hyperparameters_part1, hyperparameters_part2)

head(hyperparameters)

ile = nrow(hyperparameters)
ile_w_batchu = 250
cat(paste0("\nHow many to check: ", ile))


if(!file.exists(nazwa_konfiguracji)) {
  batch_start = 1
} else {
  tempres = data.table::fread(nazwa_konfiguracji)
  last = as.numeric(str_split(tempres$model_id, "-")[[nrow(tempres)]][1])
  if(last >= ile) { stop(paste0(last, " - Task already finished.")) }
  batch_start = as.numeric(last)+1
}

ile_batchy = ceiling(nrow(hyperparameters)/ile_w_batchu - batch_start/ile_w_batchu)

cat(paste0("\nBatch start: ", batch_start))
cat(paste0("\nHow many in batch: ", ile_w_batchu))


for (i in 1:ile_batchy) {
  batch_end = batch_start + (ile_w_batchu-1)
  if (batch_end > ile) { batch_end = ile }
  cat(paste0("\n\nProcessing batch no ", i , " of ", ile_batchy, " (", batch_start, "-", batch_end, ")"))
  
  OmicSelector_deep_learning(selected_miRNAs = selected_miRNAs, wd = getwd(), save_threshold_trainacc = 0.7, save_threshold_testacc = 0.5, hyperparameters = hyperparameters,
                             SMOTE = balanced, start = batch_start, end = batch_end, output_file = nazwa_konfiguracji, keras_threads = 30,
                             keras_epoch = 2000, keras_patience = 100, automatic_weight = F)
  
  batch_start = batch_end + 1
}

cat("[OmicSelector: TASK COMPLETED]")
sink() 
sink(type = "message")