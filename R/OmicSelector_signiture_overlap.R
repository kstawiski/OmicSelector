#' OmicSelector_signature_overlap
#'
#' A function to generate venn diagram and check the overlap between formulas.
#'
#' @param which_formulas Which formulas to check?
#' @param benchmark_csv Which benchmark to use?
#'
#' @return Object of `venn()` function which can be used for plotting venn diagram and check the overlap.
#'
#' @export
OmicSelector_signature_overlap = function(which_formulas = c("sig","cfs"), benchmark_csv = "benchmark.csv")
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
  wybrane = list()
  for (i in 1:length(which_formulas)) {
    ktora_to = match(which_formulas[i], rownames(benchmark))
    temp = as.formula(benchmark$miRy[ktora_to])
    wybrane[[rownames(benchmark)[ktora_to]]] = all.vars(temp)[-1]
  }
  require("gplots")
  temp = venn(wybrane)
  temp
}
