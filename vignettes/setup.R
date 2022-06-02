## Default repo
# Ubuntu: apt install default-jre default-jdk libmagick++-dev zlib1g-dev libglu1-mesa-dev freeglut3-dev mesa-common-dev r-cran-rgl
#r = getOption("repos")
#r["CRAN"] = "https://cran.r-project.org"
#options(repos = r)

options(repos = c(CRAN = "https://packagemanager.rstudio.com/cran/latest"))
if(Sys.info()["sysname"] == "Linux") {
  cat("This is Linux.")
  if(system("awk -F= '/^NAME/{print $2}' /etc/os-release", intern = T) == "\"Ubuntu\"") {
    cat(" Distro: UBUNTU ")
    codename = system("cat /etc/os-release | grep UBUNTU_CODENAME | cut -d = -f 2", intern = T)
    cat(codename)
    options(repos = c(CRAN = paste0("https://packagemanager.rstudio.com/cran/__linux__/",codename,"/latest")))
  }
  
}


tylko_cran = c("BiocManager","devtools","reticulate","remotes","keras","parallel")
if (length(setdiff(tylko_cran, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(tylko_cran, rownames(installed.packages())), ask = F)  }
options(Ncpus = parallel::detectCores())
#options(BioC_mirror = "https://packagemanager.rstudio.com/bioconductor")

packages = c("remotes","devtools","parallel","rlang","ps","roxygen2", "plotly", "rJava", "mice","BiocManager", "MatchIt","curl",
                       "reticulate", "kableExtra","plyr","dplyr","edgeR","epiDisplay","rsq","MASS","Biocomb","caret","dplyr",
                       "pROC","ggplot2", "doParallel", "Boruta", "spFSR", "varSelRF", "stringr", "psych", "C50", "randomForest", "doParallel",
                       "foreach","data.table", "ROSE", "deepnet", "gridExtra", "stargazer","gplots","My.stepwise","snow", "sva", "Biobase",
                       "calibrate", "ggrepel", "networkD3", "VennDiagram","RSNNS", "kernlab", "car", "PairedData",
                       "profileR","classInt","kernlab","xgboost", "keras", "tidyverse", "cutpointr","tibble","tidyr",
                       "rpart", "party", "mgcv", "GDCRNATools", "rJava", "cutpointr", "HTqPCR", "nondetects",
                       "imputeMissings", "visdat", "naniar", "stringr", "R.utils", "TCGAbiolinks", "GDCRNATools",
                       "kableExtra", "VIM", "mice", "MatchIt", "XML", "rmarkdown", "xtable", "ComplexHeatmap","circlize", "hash","RANN",
                       "BiocStyle","magick", "BiocCheck","cluster","tidyselect","ellipsis","funModeling", "mnormt","xlsx","klaR","glmnet","summarytools","shinyfullscreen","rhandsontable","shinyjs","rms","Rfast","RcppGSL","gtsummary","tableone","readr")

if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  BiocManager::install(setdiff(packages, rownames(installed.packages())), ask = F)  }

library(devtools)
library(remotes)
# Paczki z githuba
if("DMwR" %in% rownames(installed.packages()) == FALSE) { remotes::install_github("cran/DMwR", upgrade = "never") }
if("kerasformula" %in% rownames(installed.packages()) == FALSE) { remotes::install_github("rdrr1990/kerasformula", upgrade = "never") }
if("bounceR" %in% rownames(installed.packages()) == FALSE) { remotes::install_github("STATWORX/bounceR", upgrade = "never") }
if("ggbiplot" %in% rownames(installed.packages()) == FALSE) { remotes::install_github("vqv/ggbiplot", upgrade = "never") }
if("mnormt" %in% rownames(installed.packages()) == FALSE) { remotes::install_github("cran/mnormt", upgrade = "never") }
if("parsetools" %in% rownames(installed.packages()) == FALSE) { remotes::install_github("RDocTaskForce/parsetools", upgrade = "never") } 
if("testextra" %in% rownames(installed.packages()) == FALSE) { remotes::install_github("RDocTaskForce/testextra", upgrade = "never") } 
if("purrrogress" %in% rownames(installed.packages()) == FALSE) {  remotes::install_github("halpo/purrrogress", upgrade = "never") } 
if("feseR" %in% rownames(installed.packages()) == FALSE) { remotes::install_github("enriquea/feseR", upgrade = "never") }
if("autokeras" %in% rownames(installed.packages()) == FALSE) { remotes::install_github("jcrodriguez1989/autokeras", upgrade = "never") }
if("waiter" %in% rownames(installed.packages()) == FALSE) { remotes::install_github("JohnCoene/waiter", upgrade = "never") }
if("shinyjqui" %in% rownames(installed.packages()) == FALSE) { remotes::install_github("Yang-Tang/shinyjqui", upgrade = "never") }
if("RSQLite" %in% rownames(installed.packages()) == FALSE) { remotes::install_github("r-dbi/RSQLite", upgrade = "never") }

# Base dependencies
if(!'arules' %in% installed.packages()) install.packages('arules')
if(!'bit64' %in% installed.packages()) install.packages('bit64')
if(!'data.table' %in% installed.packages()) install.packages('data.table')
if(!'doParallel' %in% installed.packages()) install.packages('doParallel')
if(!'foreach' %in% installed.packages()) install.packages('foreach')
if(!'lubridate' %in% installed.packages()) install.packages('lubridate')
if(!'timeDate' %in% installed.packages()) install.packages('timeDate')

# RemixAutoML
# devtools::install_github('AdrianAntico/RemixAutoML', upgrade = FALSE, dependencies = FALSE, force = TRUE)

# tryCatch(
#         {
#             if(grepl("64", Sys.info()[["machine"]], fixed = TRUE)) {
#             # Keras
#             library(keras)
#             if(!keras::is_keras_available()) {

#               install_keras() }
#             } else { message("\n\n!!!!! If you are not running 64-bit based machine you might experience problems with keras and tensorflow that are unrelated to this package. !!!!!\n\n") }

#         },
#         error=function(cond) {
#             message(cond)
#             message("Unable to verify the correctness of keras installation. Please run keras::install_keras() later.")
#         },
#         warning=function(cond) {
#             message(cond)
#             message("Unable to verify the correctness of keras installation. Please run keras::install_keras() later.")
#         },
#         finally={

#         }
#     )

if(grepl("64", Sys.info()[["machine"]], fixed = TRUE) && !keras::is_keras_available()) { message("Keras is not installed. Please run keras::install_keras() later.") }

# OmicSelector
if("OmicSelector" %in% rownames(installed.packages()) == FALSE) { remotes::install_github("kstawiski/OmicSelector", upgrade = "never") }
message("OK! OmicSelector is installed correctly!")
