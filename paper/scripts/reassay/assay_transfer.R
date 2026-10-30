suppressPackageStartupMessages(library(pkgload))
ROOT<-normalizePath(Sys.getenv('OMICSELECTOR_PAPER_ROOT'));setwd(ROOT)
pkgload::load_all('runtime/packages/2.6.9/OmicSelector',quiet=TRUE)
source('analyses/robustness/scripts/assay_transfer.R')
