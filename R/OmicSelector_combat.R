#' OmicSelector_combat
#'
#' Use combat to fight batch effect.
#'
#' @param danex Matrix with miRNA normalized expression values (e.g. TPM, deltaCt) with miRNAs in columns and cases in rows.
#' @param metadane Metadata with `Class` and `Batch` variable.
#' @param model Model of correction (application of covariates).

#'
#' @return Batch-corrected dataset.
#'
#' @export
OmicSelector_combat = function(danex, metadane = metadane, model = c("~ Batch", "~ Batch + Class")) {

  suppressMessages(library(dplyr))
  suppressMessages(library(devtools))
  suppressMessages(library(stringr))
  suppressMessages(library(data.table))
  suppressMessages(library(tidyverse))
  suppressMessages(library(sva))
  suppressMessages(library(Biobase))
  suppressMessages(library(devtools))

  #data check
  if(table(colnames(metadane))["Class"] != 1 || length(unique(metadane$Class)) != 2) {
    stop("Metadata dataframe must contain exactly one binary 'Class' variable")
  }

  if(table(colnames(metadane))["Batch"] != 1) {
    stop("Metadata dataframe must contain nominal 'Batch' variable")
  }

  danex = as.matrix(danex)

eset <- new("ExpressionSet", exprs = t(danex))

# Covs:
pData(eset)$Batch = metadane$Batch
pData(eset)$Class = metadane$Class

# Z brakami
edata = exprs(eset)
batch = as.factor(metadane$Batch)
#batch = pheno$batch
if(model == "~ Batch") { modcombat = model.matrix(~ 1, data=metadane) }
if(model == "~ Batch + Class") { modcombat = model.matrix(~ Class, data=metadane) }

combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
rownames(combat_edata) = rownames(t(danex))

return(cbind(metadane, t(combat_edata)))
}
