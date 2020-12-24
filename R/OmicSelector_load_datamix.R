#' OmicSelector_load_datamix
#'
#' This function loads the data created in preparation phase.
#' It requires the output constructed by `OmicSelector_prepare_split` function to be placed in working directory (`wd`), thus files `mixed_train.csv`, `mixed_test.csv` and `mixed_valid.csv` have to exist in the directory.
#' For imbalanced data, the fuction can perform balancing using:
#' 1. ROSE: https://journal.r-project.org/archive/2014/RJ-2014-008/RJ-2014-008.pdf - by default we generate 10 * number of cases in orginal dataset.
#' 2. SMOTE (default): https://arxiv.org/abs/1106.1813 - by defult we use `perc.under=100` and `k=10`.
#'
#' @param wd Working directory with files for the loading.
#' @param use_smote_not_rose Set TRUE for SMOTE instead of ROSE.
#' @param smote_over Oversampling of minority class in SMOTE function (deterimes the number of cases in final dataset). See `perc.over` in `DMwR::SMOTE()`` function.
#' @param replace_smote For some analyses we may want to replace imbalanced train dataset with balanced dataset. This saved coding time in some functions.
#' @param selected_miRNAs If null - take all features staring with "hsa", if set - vector of feature names to be selected.
#'
#' @return The list of objects in the following order: train, test, valid, train_smoted, trainx, trainx_smoted, merged. (trainx contains only the miRNA data without metadata)
#'
#' @export
OmicSelector_load_datamix = function(wd = getwd(), smote_over = 10000, use_smote_not_rose = T, replace_smote = F, selected_miRNAs = NULL) {
  suppressMessages(library(plyr))
  suppressMessages(library(dplyr))
  suppressMessages(library(edgeR))
  suppressMessages(library(epiDisplay))
  suppressMessages(library(rsq))
  suppressMessages(library(MASS))
  suppressMessages(library(Biocomb))
  suppressMessages(library(caret))
  suppressMessages(library(dplyr))
  suppressMessages(library(epiDisplay))
  suppressMessages(library(pROC))
  suppressMessages(library(ggplot2))
  suppressMessages(library(DMwR))
  suppressMessages(library(ROSE))
  suppressMessages(library(gridExtra))
  suppressMessages(library(gplots))
  suppressMessages(library(devtools))
  suppressMessages(library(stringr))
  suppressMessages(library(data.table))
  suppressMessages(library(tidyverse))
  oldwd = getwd()
  setwd(wd)
  
  if(is.null(selected_miRNAs)) {
  train = dplyr::select(read.csv("mixed_train.csv", stringsAsFactors = F), starts_with("hsa"), Class)
  test = dplyr::select(read.csv("mixed_test.csv", stringsAsFactors = F), starts_with("hsa"), Class)
  valid = dplyr::select(read.csv("mixed_valid.csv", stringsAsFactors = F), starts_with("hsa"), Class) } else {
    temp = c(selected_miRNAs, "Class")
    train = dplyr::select(read.csv("mixed_train.csv", stringsAsFactors = F), temp)
    test = dplyr::select(read.csv("mixed_test.csv", stringsAsFactors = F), temp)
    valid = dplyr::select(read.csv("mixed_valid.csv", stringsAsFactors = F), temp, Class)
  }
  
  train$Class = factor(train$Class, levels = c("Control","Cancer"))
  test$Class = factor(test$Class, levels = c("Control","Cancer"))
  valid$Class = factor(valid$Class, levels = c("Control","Cancer"))

  # Wywalamy miRy z zerowa wariancja (nie musimy bo robi to poprzednia funkcja)
  #temp = train %>% dplyr::filter(Class == "Cancer")
  #temp2 = as.numeric(which(apply(temp, 2, var) == 0))
  #temp = train %>% dplyr::filter(Class == "Control")
  #temp3 = as.numeric(which(apply(temp, 2, var) == 0))
  #temp4 = unique(c(temp2, temp3))
  #if (length(temp4) > 0) {
  #  train = train[,-temp4]
  #  test = test[,-temp4]
  #  valid = valid[,-temp4]
  #}

  #train = as_tibble(train)
  #test = as_tibble(test)
  #valid = as_tibble(valid)

if(file.exists("mixed_train_balanced.csv")) {
  train_smoted = data.table::fread("mixed_train_balanced.csv")
} else {
  cat("Balanced dataset will be save as mixed_train_balanced.csv")
  if(use_smote_not_rose) {
    train_smoted = DMwR::SMOTE(Class ~ ., data = train, perc.over = smote_over,perc.under=100, k=10)
    train_smoted$Class = factor(train_smoted$Class, levels = c("Control","Cancer"))
    train_smoted = train_smoted[complete.cases(train_smoted), ]
    data.table::fwrite(train_smoted, "mixed_train_balanced.csv")
  } else {
    rosed = ROSE(Class ~ ., data = train, N = nrow(train)*10, seed = 1)
    train_smoted = rosed[["data"]]
    train_smoted$Class = factor(train_smoted$Class, levels = c("Control","Cancer"))
    train_smoted = train_smoted[complete.cases(train_smoted), ]
    data.table::fwrite(train_smoted, "mixed_train_balanced.csv")
  }
}



  if (replace_smote == T) { train = train_smoted }

  if(is.null(selected_miRNAs)) {
      trainx = dplyr::select(train, starts_with("hsa"))
      trainx_smoted = dplyr::select(train_smoted, starts_with("hsa"))
  } else {
      trainx = dplyr::select(train, selected_miRNAs)
      trainx_smoted = dplyr::select(train_smoted, selected_miRNAs)    
  }

  if(!file.exists("merged.csv")) {
    cat("Writing merged.csv for reference.")
    train2 = data.table::fread("mixed_train.csv")
    train2$mix = "train"
    test2 = data.table::fread("mixed_test.csv")
    test2$mix = "test"
    valid2 = data.table::fread("mixed_valid.csv")
    valid2$mix = "valid"
    train_balanced2 = data.table::fread("mixed_train_balanced.csv")
    train_balanced2$mix = "train_balanced"
    merged_init = rbind.fill(train2, test2, valid2, train_balanced2)
    data.table::fwrite(merged_init, "merged.csv")
  } else { merged_init = data.table::fread("merged.csv") }
  
  setwd(oldwd)
  return(list(train, test, valid, train_smoted, trainx, trainx_smoted, merged_init))
}
