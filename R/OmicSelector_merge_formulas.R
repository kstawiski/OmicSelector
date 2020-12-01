#' OmicSelector_merge_formulas
#'
#' Merge and filter formulas.
#' This function can be used to merge the formulas*.RDS created by different runs of `OmicSelector_OmicSelector()` and filter them to keep the maximum number of miRNAs.
#' This may be useful in planning qPCR validation.
#' The result of this function is `featureselection_formulas_final.RDS` file, which can be futher supplied to `OmicSelector_benchmark()`.
#'
#' @param wd Working directory with formulas*.RDS files.
#' @param max_miRNAs Maximum number of miRNAs to be selected in formulas.
#' @param add List of additional sets that should be added to the formulas.
#'
#' @return Final formulas object, also save as `featureselection_formulas_final.RDS` in working directory.
#'
#' @export
OmicSelector_merge_formulas = function(wd = getwd(), max_miRNAs = 11, add = list()) {
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
  formulas_files = list.files("temp","formulas*", all.files = T, full.names = T)
  formulas = list()
  formulas_names = character()
  for (i in 1:length(formulas_files)) {
    temp = readRDS(formulas_files[i])
    tempn = names(temp)
    formulas = c(formulas, temp)
    formulas_names = c(formulas_names, tempn)
  }

  temp = data.frame(name = formulas_names, formula = unlist(as.character(formulas)), stringsAsFactors = F)
  final = temp %>% dplyr::distinct()

  final$formula
  formulas = list()
  all = as.list(final$formula)
  names(all) =  make.names(final$name, unique = T)
  saveRDS(all, "featureselection_formulas_all.RDS")
  for(i in 1:nrow(final)){
    final$ile_miRNA[i] = length(all.vars(as.formula(final$formula[i])))-1
  }

  if(length(add) > 0) {
  for (i in 1:length(add)) {
    tempdupa = data.frame(name = names(add)[i], formula = paste0("Class ~ ", as.character(OmicSelector_create_formula(add[[i]]))[3]), ile_miRNA = 0, stringsAsFactors = F)
    final = rbind(tempdupa, final)
  }
  }

  finalold = final
  final = final %>% filter(ile_miRNA <= max_miRNAs)
  if (("fcsig" %in% final$name) == FALSE) { final = rbind(finalold[which(finalold$name=="fcsig"),], final)
  }
  if (("cfs_sig" %in% final$name) == FALSE) { final = rbind(finalold[which(finalold$name=="cfs_sig"),], final) }
  formulas_final = as.list(final$formula)
  names(formulas_final) = make.names(final$name, unique = T)
  setwd(oldwd)

  # Deduplicate before save
  f2 = unique(formulas_final)
  f2 = as.list(f2)
  for (i in 1:length(f2))
  {
        names(f2)[i] = names(formulas_final)[match(as.character(f2[[i]]),as.character(formulas_final))] # only one is matched
  }



  saveRDS(f2, "featureselection_formulas_final.RDS")
  #sink()
  #dev.off()
  fwrite(finalold, "featureselection_formulas_all.csv")
  fwrite(final, "featureselection_formulas_final.csv")
  return(formulas_final)
}
