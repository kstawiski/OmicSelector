#' OmicSelector_process_tissue_miRNA_TCGA
#'
#' Process the data downloaded from TCGA.
#'
#' @param data_folder Directory where TCGA data were downloaded.
#' @param remove_miRNAs_with_null_var Wheter to remove the miRNAs without any expression? Default: True
#'
#' @export
OmicSelector_process_tissue_miRNA_TCGA = function(data_folder = getwd(), remove_miRNAs_with_null_var = T) {
  suppressMessages(library(plyr))
  suppressMessages(library(data.table))
  suppressMessages(library(dplyr))
  suppressMessages(library(edgeR))
  suppressMessages(library(naniar))
  suppressMessages(library(visdat))
  suppressMessages(library(stringr))
  suppressMessages(library(foreach))
  suppressMessages(library(doParallel))
  suppressMessages(library(imputeMissings))
  pliki = list.files(data_folder, "*.csv", all.files = T, full.names = T)[startsWith(list.files(data_folder, "*.csv"), "miRNA_")]


  counts = fread(pliki[1])
  counts = OmicSelector_correct_miRNA_names(counts)
  for (i in 2:length(pliki)) {
    temp = fread(pliki[i])
    temp = OmicSelector_correct_miRNA_names(temp)
    counts = rbind.fill(counts, temp)
  }

  str(counts)

  na_count <-sapply(counts, function(y) sum(length(which(is.na(y)))))
  na_count <- data.frame(na_count)
  na_count$miR = row.names(na_count)

  projects = read.csv(paste0(data_folder,"/projects.csv"))
  counts$primary_site = projects$primary_site[match(counts$project_id, projects$project_id)]
  counts$disease_type = projects$disease_type[match(counts$project_id, projects$project_id)]
  counts$site = projects$name[match(counts$project_id, projects$project_id)]

  dane_counts = dplyr::select(counts, starts_with("hsa-"))
  #colnames(dane_counts) = make.names(colnames(dane_counts),unique = T)

  na_count <-sapply(dane_counts, function(y) sum(length(which(is.na(y)))))
  na_count <- data.frame(na_count)
  na_count$miR = row.names(na_count)

  head(na_count[order(-na_count$na_count),])

  # usuwamy miRNA z więcej niż 10% missing
  dane_counts = dane_counts[ , colSums(is.na(dane_counts)) < (0.1*nrow(dane_counts))]

  # usuwamy miRNA z jakimikolwiek missing
  #dane_counts = dane_counts[ , colSums(is.na(dane_counts)) == 0]

  # imputujemy resztę
  #suppressMessages(library(mice))
  #temp1 = parlmice(dane_counts, m=1, n.core = detectCores()-1, seed = 1, n.imp.core = 50, cl.type = "FORK")
  #temp2 = temp1$data
  #dane_counts_missing = dane_counts
  #dane_counts = complete(temp1)
  # trwało wieki.. imputowanie średnią
  dane_counts = impute(dane_counts, method = "median/mode", flag = FALSE)

  na_count <-sapply(dane_counts, function(y) sum(length(which(is.na(y)))))
  na_count <- data.frame(na_count)
  na_count$miR = row.names(na_count)

  head(na_count[order(-na_count$na_count),])

  dane_p = dplyr::select(counts, -starts_with("hsa-"))

  fwrite(cbind(dane_p, dane_counts), "tissue_miRNA_counts.csv")


  dane_counts = t(dane_counts)
  colnames(dane_counts) = counts$sample

  # Normalizacja do TPM
  dane3 = DGEList(counts=dane_counts, genes=data.frame(miR = rownames(dane_counts)), samples = dane_p)
  tpm = cpm(dane3, normalized.lib.sizes=F, log = F, prior.count = 0.001)
  tpm = tpm + 0.001
  tpm = log10(tpm)
  ttpm = t(tpm)

  # Usuwamy miRNA z 0 warancją
  if(remove_miRNAs_with_null_var) {
    ttpm <- ttpm[, - as.numeric(which(apply(ttpm, 2, var) == 0))] }

  fwrite(cbind(dane_p, ttpm), "tissue_miRNA_logtpm.csv")




}
