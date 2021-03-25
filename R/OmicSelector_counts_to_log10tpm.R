#' OmicSelector_counts_to_log10tpm
#'
#' Counts to log-transformed TPM-normalized counts.
#' The funcction support additional filter, i.e. it can leave miRNAs that appear in minimum of `filtr_minimalcounts` counts in `filtr_howmany` samples.
#' Usage of the filter like that can assure that the miRNAs selected will be detectable in qPCR.
#'
#' @param danex Matrix with miRNA counts with miRNAs in columns and cases in rows.
#' @param metadane Metadata with `Class` variable.
#' @param  ids Unique identifier of samples.
#' @param filtr If expression filter should be used.
#' @param filtr_minimalcounts How many counts?
#' @param filtr_howmany In how many samples? (Please provide a percentage or proportion, e.g. 1/2 or 1/3).
#' @param increment Increment added to TPM values before log-transformation (usually: 0.001, so -3 will be equalt to lack of expression).
#'
#' @return Normalized counts in `ttpm` format. Please note, that the function also saves `TPM_DGEList_filtered.rds` to working directory, which is a DGEList object that can be used in packages like edgeR.
#'
#' @export

OmicSelector_counts_to_log10tpm = function(danex, metadane = metadane, ids = metadane$ID, filtr = T,
                                 filtr_minimalcounts = 10, filtr_howmany = 1/2, increment = 0.001) {
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

  #data check

  #check if all columns are numerical and with names starting as hsa
  for(i in colnames(danex)) {
    if(!is.numeric(danex[, i])) {
      stop("Please provide a dataframe with only numeric variables")
    }
    if(!startsWith(i, "hsa")){
      stop("Please provide only microRNA expression data (column names starting with hsa)")
    }
  }
  if(table(colnames(metadane))["Class"] != 1 || length(unique(metadane$Class)) != 2) {
    stop("Metadata dataframe must contain exactly one binary 'Class' variable")
  }

  danex = as.matrix(danex)

  dane_counts = t(danex)
  colnames(dane_counts) = ids
  mode(dane_counts) = 'numeric'

  # Normalizacja do TPM
  dane3 = DGEList(counts=dane_counts, genes=data.frame(miR = rownames(dane_counts)), samples = metadane)
  tpm = cpm(dane3, normalized.lib.sizes=F, log = F, prior.count = 0.001)
  tpm = tpm + increment
  tpm = log10(tpm)
  ttpm = t(tpm)
  saveRDS(dane3,"TPM_DGEList.rds")
  cat("\nDGEList unfiltered object with TPM was saved as TPM_DGEList.rds.")

  if (filtr == F) {
    cat("\nReturned data are log10(TPM).")
    return(ttpm)
  } else {
    zostaw = vector()
    zostaw[1] = TRUE
    for (i in 1:nrow(dane3))
    {
      temp = as.numeric(dane3$counts[i,])
      if (sum(temp < filtr_minimalcounts) > (filtr_howmany)*length(temp)) { zostaw[i] = FALSE } else { zostaw[i] = TRUE }
    }

    dane4 = dane3[zostaw,]
    ttpm_features = ttpm[,which(zostaw)]

    cat("\nDGEList filtered object with TPM was saved as TPM_DGEList_filtered.rds.")
    cat(paste0("\n(After filtering) miRNAs left: ", sum(zostaw), " | filtered out: ", sum(zostaw==FALSE), "."))
    saveRDS(dane4,"TPM_DGEList_filtered.rds")
    cat("\nReturned data are log10(TPM).")
    return(ttpm_features)
  }

}
