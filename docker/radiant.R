#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if("radiant" %in% rownames(installed.packages()) == FALSE) {
  devtools::install_github("radiant-rstats/radiant")
}

if(length(args) > 0) {
if(!dir.exists("/radiant-data")) { dir.create("/radiant-data") }
unlink("/radiant-data/*")
options(radiant.init.data = list.files(path = "/radiant-data", full.names = TRUE))
setwd(paste0("/OmicSelector/",args[1],"/"))
pliki = list.files(path = paste0("/OmicSelector/",args[1],"/"), pattern="*.csv", recursive = F)
# pliki2 = list.files(path = paste0("/OmicSelector/",args[1],"/"), pattern="merged(.*).csv", recursive = T)
# pliki = c(pliki1, pliki2)
id = args[1]
for(i in 1:length(pliki)){
  try({ temp = read.csv(pliki[i])
  save(temp, file = paste0("/radiant-data/",id, "_", make.names(pliki[i]), ".rda")) })
}
  setwd("/radiant-data")
  # try({ procesy = system("lsof -i | grep 3839", intern = T)
  # pid = strsplit(procesy, " ")[[1]][strsplit(procesy, " ")[[1]] != ""][2]
  # system(paste0("kill -9 ", pid)) })
  # system('screen -dmS radiant Rscript -e "library(radiant); library(data.table); radiant_url(port = 3839);"')

} else {
  setwd("/radiant-data")
  # system('screen -dmS radiant Rscript -e "library(radiant); library(data.table); radiant_url(port = 3839);"')
  print("No analysisid.")
}