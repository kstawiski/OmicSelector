#' OmicSelector_prepare_split
#'
#' Create split required for `OmicSelector_OmicSelector()` and all the following functions. It is obligatory to use it.
#' The function devides the dataset into training, testing and validation set. Be default (as `train_proc=0.6`) 60 perc. of cases will be assigned to trainining datset.
#' The rest is devided into testing and validation dataset in half, ending in 60 perc. of cases in training dataset, 20 perc. of cases in testing dataset and 20 perc. of cases in validation dataset.
#' Metadata have to have `Class` variable, with `Cancer` and `Control` values.
#'
#' @param metadane Metadata of cases. Must contain `Class` variable with `Cancer` and `Control` values.
#' @param ttpm Normalized counts used (primary data for the rest of the analysis).
#' @param train_proc What perc. should be kept in training dataset?
#'
#' @return The mixed dataset is return. In working directory mixed_train.csv, mixed_test.csv and mixed_valid.csv are saved. This is a crucial step in data preprocessing.
#'
#' @export
OmicSelector_prepare_split = function(metadane = metadane, ttpm = ttpm_features, train_proc = 0.6)
{
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
  tempp = cbind(metadane, ttpm)
  # Podzial - http://topepo.github.io/caret/data-splitting.html#simple-splitting-based-on-the-outcome

  #data check

  #check if all columns are numerical and with names starting as hsa
  for(i in colnames(ttpm)) {
    if(!is.numeric(ttpm[, i])) {
      stop("Please provide a dataframe with only numeric variables")
    }
    if(!startsWith(i, "hsa")){
      stop("Please provide only microRNA expression data (column names starting with hsa)")
    }
  }

  metadane <- data.frame(metadane)

  if(table(colnames(metadane))["Class"] != 1 || length(unique(metadane$Class)) != 2) {
    stop("Metadata dataframe must contain exactly one binary 'Class' variable")
  }


  set.seed(1)
  mix = rep("unasign", nrow(tempp))
  suppressMessages(library(caret))
  train.index <- createDataPartition(tempp$Class, p = train_proc, list = FALSE)
  mix[train.index] = "train"
  train = tempp[train.index,]


  tempp2 =  tempp[-train.index,]
  train.index <- createDataPartition(tempp2$Class, p = .5, list = FALSE)
  test = tempp2[train.index,]
  valid = tempp2[-train.index,]
  write.csv(train, "mixed_train.csv",row.names = F)
  write.csv(test, "mixed_test.csv",row.names = F)
  write.csv(valid, "mixed_valid.csv",row.names = F)

  train$mix = "train"
  test$mix = "test"
  valid$mix = "valid"

  mixed = rbind(train,test,valid)

  metadane2 = dplyr::select(mixed, -starts_with("hsa"))
  ttpm2 = dplyr::select(mixed, starts_with("hsa"))

  mixed2 = cbind(metadane2, ttpm2)
  write.csv(mixed2, "mixed.csv",row.names = F)
  cat("\nSaved 3 sets as csv in working directory. Retruned mixed dataset.")
  return(mixed)
}
