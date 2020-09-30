#' OmicSelector_get_features_from_benchmark
#'
#' Get which features are used by which method in the benchmark.
#'
#' @param benchmark_csv Path to benchmark csv.
#' @param method Method of interest.
#'
#' @return Vector of miRNAs.
#'
#' @export
OmicSelector_get_features_from_benchmark = function(benchmark_csv = "benchmark1578990441.6531.csv", method = "fcsig")
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
  benchmark = read.csv(benchmark_csv, stringsAsFactors = F)
  rownames(benchmark) = make.names(benchmark$method, unique = T)
  return(all.vars(as.formula(benchmark$miRy[which(rownames(benchmark) == method)]))[-1])
}
