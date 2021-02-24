# Parameters:
balanced = F
if(file.exists("var_deeplearning_balanced.txt")) { balanced = as.logical(readLines("var_deeplearning_balanced.txt", warn = F)) }
autoencoders = F
if(file.exists("var_deeplearning_autoencoders.txt")) { autoencoders = as.numeric(readLines("var_deeplearning_autoencoders.txt", warn = F)) }
keras_threads = 30
if(file.exists("var_deeplearning_keras_threads.txt")) { keras_threads = as.numeric(readLines("var_deeplearning_keras_threads.txt", warn = F)) }
if(file.exists("var_deeplearning_selected.txt")) { selected_miRNAs = readLines("var_deeplearning_selected.txt", warn = F) }
if(selected_miRNAs != "all")
{
  library(OmicSelector)
  input_formulas = readRDS("featureselection_formulas_final.RDS")
  miRNAs = all.vars(as.formula(input_formulas[[selected_miRNAs]]))[-1];
  if(length(miRNAs)>0) { selected_miRNAs = miRNAs }
} else { selected_miRNAs = colnames(data.table::fread("mixed_train.csv"))[startsWith(colnames(data.table::fread("mixed_train.csv")),"hsa")] }


# Code:
nazwa_konfiguracji = "deeplearning.csv"
options(warn = -1)
if(file.exists("task.log")) { file.remove("task.log") }
con <- file("task.log")
sink(con, append=TRUE)
sink(con, append=TRUE, type = "message")
library(OmicSelector); set.seed(1);
if(!dir.exists("/OmicSelector/temp")) { dir.create("/OmicSelector/temp") }
OmicSelector_load_extension("deeplearning")
library(data.table)

# Load check
current = 1
while(current > 0.5) { 
  load = strsplit(system("cat /proc/loadavg", intern = T)," ")
  max = parallel::detectCores()
  current = as.numeric(load[[1]][1])/max
  if(current > 0.5) {
    try({ OmicSelector_log(paste0("Current server load: ", round(current*100,2), "% exceeds the threshold of 50%. The job waiting for resources to start...\n"),"task.log") })
    Sys.sleep(15);
  } else { cat(paste0("Current server load: ", round(current*100,2), "%. The job is starting...\n")); }}

try({
gpu = "none"
gpu = system("lspci -vnnn | perl -lne 'print if /^\\d+\\:.+(\\[\\S+\\:\\S+\\])/' | grep VGA", intern = T) 
if(grepl("NVIDIA", gpu)) {
gpu_util = 100
gpu_mem = 1
gpu_memu = 1

while(gpu_util > 50 && (gpu_memu/gpu_mem) > 0.5) {

  try({
  gpu_util = as.numeric(system("nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader,nounits", intern = T))
  gpu_mem = as.numeric(system("nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits", intern = T))
  gpu_memu = as.numeric(system("nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits", intern = T))
  try({ OmicSelector_log(paste0("\nGPU | Util: ", gpu_util, "% | Memory: ", round(gpu_memu/gpu_mem, 4)*100, "%"),"task.log") })
  })
  
  try({ OmicSelector_log(paste0("\nWaiting for resources to start the task..."),"task.log") })
  Sys.sleep(10)
}
try({ OmicSelector_log(paste0("\nTask is starting..."),"task.log") })
}})

# Data
dane = OmicSelector_load_datamix()
if(balanced == F) {
  t = data.table::fread("mixed_train.csv")
} else {
  t = data.table::fread("mixed_train_balanced.csv")
}


# autoencoder ma softmax na deep feature
if(autoencoders == 1) { hyperparameters_part1 = expand.grid(layer1 = seq(2,10, by = 1), layer2 = c(0), layer3 = c(0),
                                    activation_function_layer1 = c("relu","sigmoid","selu"), activation_function_layer2 = c("relu"), activation_function_layer3 = c("relu"),
                                    dropout_layer1 = c(0, 0.1), dropout_layer2 = c(0), dropout_layer3 = c(0),
                                    layer1_regularizer = c(T,F), layer2_regularizer = c(F), layer3_regularizer = c(F),
                                    optimizer = c("adam","rmsprop","sgd"), autoencoder = c(0,-7,7), balanced = balanced, formula = as.character(OmicSelector_create_formula(selected_miRNAs))[3], scaled = c(T,F),
                                    stringsAsFactors = F)
hyperparameters_part2 = expand.grid(layer1 = seq(3,11, by = 2), layer2 = c(seq(3,11, by = 2)), layer3 = c(seq(0,11, by = 2)),
                                    activation_function_layer1 = c("relu","sigmoid","selu"), activation_function_layer2 = c("relu","sigmoid","selu"), activation_function_layer3 = c("relu","sigmoid","selu"),
                                    dropout_layer1 = c(0, 0.1), dropout_layer2 = c(0), dropout_layer3 = c(0),
                                    layer1_regularizer = c(T,F), layer2_regularizer = c(F), layer3_regularizer = c(F),
                                    optimizer = c("adam","rmsprop","sgd"), autoencoder = c(0,-7,7), balanced = balanced, formula = as.character(OmicSelector_create_formula(selected_miRNAs))[3], scaled = c(T,F),
                                    stringsAsFactors = F)
hyperparameters = rbind(hyperparameters_part1, hyperparameters_part2) 
} else { 
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
hyperparameters = rbind(hyperparameters_part1, hyperparameters_part2) }

# if quick scan
if(autoencoders == 2) {
  hyperparameters = expand.grid(layer1 = seq(2,10, by = 1), layer2 = c(0), layer3 = c(0),
                                    activation_function_layer1 = c("relu","sigmoid","selu"), activation_function_layer2 = c("relu"), activation_function_layer3 = c("relu"),
                                    dropout_layer1 = c(0, 0.1), dropout_layer2 = c(0), dropout_layer3 = c(0),
                                    layer1_regularizer = c(T,F), layer2_regularizer = c(F), layer3_regularizer = c(F),
                                    optimizer = c("adam","rmsprop","sgd"), autoencoder = c(0,-7,7), balanced = balanced, formula = as.character(OmicSelector_create_formula(selected_miRNAs))[3], scaled = c(T,F),
                                    stringsAsFactors = F)
}

# head(hyperparameters)

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
  try({ OmicSelector_log(paste0("\n\nProcessing batch no ", i , " of ", ile_batchy, " (", batch_start, "-", batch_end, ")"),"task.log") })

  OmicSelector_deep_learning(selected_miRNAs = selected_miRNAs, wd = getwd(), save_threshold_trainacc = 0.7, save_threshold_testacc = 0.5, hyperparameters = hyperparameters,
                             SMOTE = balanced, start = batch_start, end = batch_end, output_file = nazwa_konfiguracji, keras_threads = 30,
                             keras_epoch = 2000, keras_patience = 100, automatic_weight = F)
  
  batch_start = batch_end + 1
}

# Merge:
lista_plikow = list.files(".", pattern = "^deeplearning.*.csv$")
library(plyr)
wyniki = data.frame()
for(i in 1:length(lista_plikow)) { temp = data.table::fread(lista_plikow[i]); wyniki = rbind.fill(wyniki, temp); }
data.table::fwrite(wyniki, "merged_deeplearning.csv")

wyniki$metaindex = (wyniki$training_Accuracy + wyniki$test_Accuracy + wyniki$valid_Accuracy) / 3
wyniki$metaindex2 = (wyniki$test_Accuracy + wyniki$valid_Accuracy) / 2
library(dplyr)
wynikitop = wyniki %>% arrange(desc(metaindex)) %>% filter(worth_saving == TRUE)
if(nrow(wynikitop)<1000) { max_top = nrow(wynikitop) } else { max_top = 1000}
wynikitop = wynikitop[1:max_top,]
data.table::fwrite(wynikitop, "merged_deeplearning_top.csv")
data.table::fwrite(as.data.frame(wynikitop$name), "merged_deeplearning_names.csv")

cat("[OmicSelector: TASK COMPLETED]")
sink() 
sink(type = "message")
try({ OmicSelector_log("[OmicSelector: TASK COMPLETED]","task.log") })