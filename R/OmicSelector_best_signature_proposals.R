#' OmicSelector_best_signature_proposals
#'
#' Propose the best signture based on benchamrk methods.
#' This function calculated the `metaindex` value which is the harmonic mean of accuracy on train, test and validation dataset.
#' In the next step, it sorts the miRNA sets based on `metaIndex1` score. The first row in resulting data frame is the winner miRNA set.
#'
#' @param benchmark_csv Path to benchmark csv.
#' @param without_train One can argue, that accuracy on training dataset should not be used in calculation of metaIndex. By setting this to TRUE, you can calculate it without train dataset.
#'
#' @return The benchmark sorted by metaIndex. First row is the best performing miRNA set.
#'
#' @export
#'
OmicSelector_best_signature_proposals = function(benchmark_csv = "benchmark1578929876.21765.csv", without_train = F){
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
  benchmark = read.csv(benchmark_csv, stringsAsFactors = F)
  rownames(benchmark) = make.names(benchmark$method, unique = T)
  
  trainacc = dplyr::select(benchmark, ends_with("_train_Accuracy"))
  trainacc_metaindex = rowMeans(trainacc)

  testacc = dplyr::select(benchmark, ends_with("_test_Accuracy"))
  testacc_metaindex = rowMeans(testacc)

  validacc = dplyr::select(benchmark, ends_with("_valid_Accuracy"))
  validacc_metaindex = rowMeans(validacc)

  acc = data.frame(trainacc_metaindex, testacc_metaindex, validacc_metaindex)
  temp = t(acc)
  acc$metaindex = psych::harmonic.mean(temp)
  
  if (without_train == T) { 
    acc = data.frame(testacc_metaindex, validacc_metaindex)
    temp = t(acc)
    acc$metaindex = psych::harmonic.mean(temp)
   }
  acc$method = rownames(benchmark)
  acc$miRy = benchmark$miRy
  rownames(acc) = make.names(benchmark$method, unique = T)
  acc = acc %>% arrange(desc(metaindex))
  return(acc)
}
