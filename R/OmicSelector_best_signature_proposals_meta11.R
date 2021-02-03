#' OmicSelector_best_signature_proposals_meta11
#'
#' Propose the best signture based on benchamrk methods.
#' This function calculated the `metaIndex11` value which is the Youden-like score on validation set (the only one that was never used in any section of the pipeline).
#' Formula: `metaIndex11 = validation sensitivitiy + validation specificity - 1`
#' In the next step, it sorts the miRNA sets based on `metaIndex11` score. The first row in resulting data frame is the winner miRNA set.
#'
#' @param benchmark_csv Path to benchmark csv.
#'
#'
#' @return The benchmark sorted by metaIndex. First row is the best performing miRNA set.
#'
#' @export
OmicSelector_best_signature_proposals_meta11 = function(benchmark_csv = "benchmark1578929876.21765.csv"){
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
  temp1 =  dplyr::select(benchmark, ends_with("_valid_Sensitivity"))
  temp2 = dplyr::select(benchmark, ends_with("_valid_Specificity"))
  #acc = dplyr::select(benchmark, ends_with("_train_Accuracy"), ends_with("_test_Accuracy"),ends_with("_valid_Accuracy") )
  #if (without_train == T) { acc = dplyr::select(benchmark, ends_with("_test_Accuracy"),ends_with("_valid_Accuracy")) }
  temp1$temp = rowMeans(temp1)
  temp2$temp = rowMeans(temp2)
  acc = benchmark
  acc$metaindex = temp1$temp + temp2$temp - 1
  acc$method = rownames(benchmark)
  acc$miRy = benchmark$miRy
  rownames(acc) = make.names(benchmark$method, unique = T)
  acc = acc %>% arrange(desc(metaindex))
  return(acc)
}
