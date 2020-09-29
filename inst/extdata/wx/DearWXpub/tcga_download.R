
#set command path
#setwd('/current_absolute_path/')

#' Load functions
source("Module_A.R")

#' set data saving path
sPath1 <- "./TCGA_DATAS"

#cancer type
type <- c("BLCA", "BRCA", "COAD", "HNSC", "KICH", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PRAD", "THCA")

#download rna-seq data
for(x in type){
  DownloadRNASeqData(cancerType = x,
                     assayPlatform = "gene.normalized_RNAseq",
                     saveFolderName = sPath1)
}
