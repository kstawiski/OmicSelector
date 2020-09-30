#' OmicSelector_correct_miRNA_names
#'
#' Sometinmes, when using the dataset mapped to previous versions of miRbase we may get false mismatches due to changes in terminology.
#' This function uses latest version of miRbase to correct all old miRNA names to new one.
#'
#' @param temp Dataset with miRNA names in columns.
#' @param correct_dots Boolean variable to correct the names after correction of dots to hyphens. This tries to compensate the effect of make.names() function.
#'
#' @return Corrected dataset.
#'
#' @export
OmicSelector_correct_miRNA_names = function(temp, species = "hsa", correct_dots = T) {
  suppressMessages(library(foreach))
  suppressMessages(library(doParallel))
  suppressMessages(library(dplyr))
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
  miRbase_aliasy = fread("ftp://mirbase.org/pub/mirbase/CURRENT/aliases.txt.gz")
  colnames(miRbase_aliasy) = c("MIMAT","Aliasy")
  miRbase_aliasy_hsa = miRbase_aliasy %>% filter(str_detect(Aliasy, paste0(species,"*"))) %>% filter(str_detect(MIMAT, "MIMAT*"))

  #setup parallel backend to use many processors
  cores=detectCores()
  cl <- makePSOCKcluster(useXDR = TRUE, cores-1) #not to overload your computer
  suppressMessages(library(doParallel))
   registerDoParallel(cl)
  # on.exit(stopCluster(cl))

  temp2 = colnames(temp)
  final <- foreach(i=1:length(temp2), .combine=c) %dopar% {
    #for(i in 1:length(temp2)) {
    if(correct_dots) { naz = gsub("\\.", "-", temp2[i]) } else { naz = temp2[i] } # correct dots to hyphens
    suppressMessages(library(data.table))
    suppressMessages(library(stringr))
    for (ii in 1:nrow(miRbase_aliasy_hsa)) {
      temp3 = str_split(as.character(miRbase_aliasy_hsa[ii,2]), ";")
      temp4 = temp3[[1]]
      temp4 = temp4[temp4 != ""]
      if(naz %in% temp4) { naz = temp4[length(temp4)] }
    }




    naz



  }

  colnames(temp) = final
  stopCluster(cl)
  return(temp)
}
