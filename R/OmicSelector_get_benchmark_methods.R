#' OmicSelector_get_benchmark_methods
#'
#' Get methods checked in benchmark.
#'
#' @param benchmark_csv Path to benchmark csv file.
#' @return Vector of feature selection methods checked.
#'
#' @export
OmicSelector_get_benchmark_methods = function(benchmark_csv = "benchmark1578929876.21765.csv"){
  suppressMessages(library(limma))
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
  benchmark$method = rownames(benchmark)
  temp = dplyr::select(benchmark, ends_with("_valid_Accuracy"))
  metody = strsplit2(colnames(temp), "_")[,1]
  return(metody)
}
