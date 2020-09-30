#' OmicSelector_differential_expression_ttest
#'
#' The variable performes standard differential expression analysis using unpaired t-test with BH and Bonferonni correction.
#' It requires `ttpm_features` object, which is e.g. a matrix of log-transformed TPM-normalized miRNAs counts with miRNAs placed in `ttpm` design (i.e. columns and cases placed as rows).
#' Classess should be passed as `classes` and this should be a vector of length equal to number of rows in `ttpm_polfiltrze` and contain only "Cancer" or "Control" labels!!
#' The function returns the miRNAs sorted by BH-corrected p-value.
#'
#' @param ttpm_features matrix of log-transformed TPM-normalized miRNAs counts or other feature matrix with miRNAs/features placed in `ttpm` design (i.e. columns and cases placed as rows)
#' @param classes vector describing label for each case. It should contain only "Cancer" and "Control" labeles!!!!
#' @param mode use 'logtpm' for log(TPM) data or 'deltact' for qPCR deltaCt values. This parameters sets how the fold-change is calculated. Setting it to "auto" will try to read settings from var_type.txt (used in docker).
#'
#'
#' @return Data frame with results.
#'
#' @export
OmicSelector_differential_expression_ttest = function(ttpm_features, classes, mode = "auto")
{

  # obsluga auto
  if(mode == "auto") {
    if(file.exists("var_type.txt"))
    { type = readLines("var_type.txt", warn = F)
     mode = as.character(type) } else {
       mode = "logtpm" # ustaw defaultowy
     }
  }


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
  wyniki = data.frame(miR = as.character(colnames(ttpm_features)))
  ttpm_features = as.data.frame(ttpm_features)
  classes = as.character(classes) # dopilnujemy, zeby sam faktoryzowal

  suppressMessages(library(dplyr))

  for (i in 1:length(colnames(ttpm_features)))
  {
    # Filtry
    #wyniki[i,"jakieszero"] = ifelse(sum(ttpm_features[,i] == 0) > 0, "tak", "nie")

    # Średnia i SD
    wyniki[i,paste0("mean ",mode)] = mean(ttpm_features[,i])
    wyniki[i,paste0("median ",mode)] = median(ttpm_features[,i])
    wyniki[i,paste0("SD ",mode)] = sd(ttpm_features[,i])

    # Cancer
    tempx = ttpm_features[classes == "Cancer",]
    wyniki[i,"cancer mean"] = mean(tempx[,i])
    wyniki[i,"cancer median"] = median(tempx[,i])
    wyniki[i,"cancer SD"] = sd(tempx[,i])

    # Cancer
    tempx = ttpm_features[classes != "Cancer",]
    wyniki[i,"control mean"] = mean(tempx[,i])
    wyniki[i,"control median"] = median(tempx[,i])
    wyniki[i,"control SD"] = sd(tempx[,i])

    # DE
    temp = t.test(ttpm_features[,i] ~ as.factor(classes))
    if(mode == "logtpm") {
    fc = (temp$estimate[1] - temp$estimate[2])
    wyniki[i,"log10FC (subtr estim)"] = fc
    wyniki[i,"log10FC"] = wyniki[i,"cancer mean"] - wyniki[i,"control mean"]
    wyniki[i,"log2FC"] = wyniki[i,"log10FC"] / log10(2)

    revfc = (temp$estimate[2] - temp$estimate[1])
    wyniki[i,"reversed_log10FC"] = revfc
    wyniki[i,"reversed_log2FC"] = wyniki[i,"reversed_log10FC"] / log10(2)
    #wyniki[i,"log2FC"] = log2(fc)
    }

    if(mode == "deltact") {
    fc = 2 ^ (temp$estimate[1] - temp$estimate[2])
    wyniki[i,"FC"] = fc
    wyniki[i,"log10FC"] = log10(fc)
    wyniki[i,"log2FC"] = log2(fc)

    revfc = 2 ^ (temp$estimate[2] - temp$estimate[1])
    wyniki[i,"reversed_log10FC"] = log10(fc)
    wyniki[i,"reversed_log2FC"] = log2(fc)
    #wyniki[i,"log2FC"] = log2(fc)
    }

    wyniki[i,"p-value"] = temp$p.value
  }
  wyniki[,"p-value Bonferroni"] = p.adjust(wyniki$`p-value`, method = "bonferroni")
  wyniki[,"p-value Holm"] = p.adjust(wyniki$`p-value`, method = "holm")
  # wyniki[,"-log10(p-value Bonferroni)"] = -log10(p.adjust(wyniki$`p-value`, method = "bonferroni"))
  wyniki[,"p-value BH"] = p.adjust(wyniki$`p-value`, method = "BH")

  wyniki$miR = as.character(wyniki$miR)
  return(wyniki %>% arrange(`p-value BH`))
}
