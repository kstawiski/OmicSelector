#' OmicSelector_OmicSelector
#'
#' Main function of the package. The aim of this function is to perform feature selection using multiple methods and to create formulas for benchmarking.
#' It loads the data from working directory. The output is mainly created in files in working directory. Log and temporary files are placed in created `temp` subfolder.
#' This package offers about 60 feature selection methods. Which methods will be check by this function is defined by `m` parameter.
#' Pearls about the methods:
#'
#' - Sig = miRNAs with p-value <0.05 after BH correction (DE using t-test)
#' - Fcsig = sig + absolute log2FC filter (included if abs. log2FC>1)
#' - Cfs = Correlation-based Feature Selection for Machine Learning (more: https://www.cs.waikato.ac.nz/~mhall/thesis.pdf)
#' - Classloop = Classification using different classification algorithms (classifiers) with the embedded feature selection and using the different schemes for the performance validation (more: https://rdrr.io/cran/Biocomb/man/classifier.loop.html)
#' - Fcfs = CFS algorithm with forward search (https://rdrr.io/cran/Biocomb/man/select.forward.Corr.html)
#' - MDL methods = minimal description length (MDL) discretization algorithm with different a method of feature ranking or feature selection (AUC, SU, CorrSF) (more: https://rdrr.io/cran/Biocomb/man/select.process.html)
#' - bounceR = genetic algorithm with componentwise boosting (more: https://www.statworx.com/ch/blog/automated-feature-selection-using-bouncer/)
#' - RandomForestRFE = recursive feature elimination using random forest with resampling to assess the performance. (more: https://topepo.github.io/caret/recursive-feature-elimination.html#resampling-and-external-validation)
#' - GeneticAlgorithmRF (more: https://topepo.github.io/caret/feature-selection-using-genetic-algorithms.html)
#' - SimulatedAnnealing =  makes small random changes (i.e. perturbations) to an initial candidate solution (more: https://topepo.github.io/caret/feature-selection-using-simulated-annealing.html)
#' - Boruta (more: https://www.jstatsoft.org/article/view/v036i11/v36i11.pdf)
#' - spFSR = simultaneous perturbation stochastic approximation (SPSA-FSR) (more: https://arxiv.org/abs/1804.05589)
#' - varSelRF = using the out-of-bag error as minimization criterion, carry out variable elimination from random forest, by successively eliminating the least important variables (with importance as returned from random forest). (more: https://www.ncbi.nlm.nih.gov/pubmed/16398926)
#' - WxNet = a neural network-based feature selection algorithm for transcriptomic data (more: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6642261/)
#' - Step = backward stepwise method of feature selection based on logistic regression (GLM, family = binomial) using AIC criteria (stepAIC) and functions from My.stepwise package (https://cran.r-project.org/web/packages/My.stepwise/index.html)
#'
#' For more detailed defitions please see the tutorial vignette.
#'
#' @param wd Working directory with data (`mixed_train.csv`, `mixed_test.csv` and `mixed_validation.csv` as created by `OmicSelector_prepare_split` have to be present).
#' @param m Methods of feature selection to be performed. This has to be a vector of integers with minimum of 1 and maximum of 70. For the definition of numbers please see the vignette.
#' @param max_iterations Maximum number of iterations in selected methods. Setting this too high may results in very long comupting time.
#' @param code_path A folder where the python external scripts are placed (especially for WxNet method). By default the additional code is provided in the package.
#' @param register_parallel Where to use parallel processing to speed up computing time. Seting it to FALSE may aid in debuging.
#' @param clx This parameter may be used for passing the already register computing cluster (created and registered with `doParallel` tools). This may lower the computing time by saving the time to register new cluster.
#' @param stamp A character vector or timestamp used for marking the output files.
#' @param prefer_no_features Maximum number of miRNAs that can be selected by the tools if the method allows for that.
#' @param conda_path Patch to "conda" bindary used for executing python scripts.
#' @param debug Gives additional debug information (saves .rdata after feature selection is completed, prints formulas to log)
#' @param timeout_sec Timeout after the method is terminated if not finished. It may be useful to keep the long methods limited, not to wait ethernity for the results.
#' @param type Parameter 'mode' forwarded to OmicSelector_differential_expression_ttest which is essential in many feature selection methods. Note that if 'var_type.txt' file exists in the working directory it is superior to the value set directly in function - as designed for GUI. Please refer to OmicSelector_differential_expression_ttest manual to understand how 'mode' worOmicSelector_
#'
#' @return The list of selected formulas. Note that, due to purpose of this package `OmicSelector_merge_formulas` may be a better option to get the output of processes run by this function.
#' @examples
#' # NOT RUN: (to speed up check, but this is a valid example for your real time projects)
#' # suppressMessages(library(foreach))
#' # suppressMessages(library(doParallel))
#' # suppressMessages(library(parallel))
#' # suppressMessages(library(doParallel))
#' # m = 1:56 # which methods to check?
#' # cl <- makePSOCKcluster(useXDR = TRUE, 5) # 5 threds by default
#' # doParallel:: registerDoParallel(cl)
#' # iterations = length(m)
#' # pb <- txtProgressBar(max = iterations, style = 3)
#' # progress <- function(n) setTxtProgressBar(pb, n)
#' # opts <- list(progress = progress)
#' # foreach(i = m, .verbose = TRUE, .options.snow = opts) %dopar%
#' # {
#' # suppressMessages(library(OmicSelector))
#' # setwd("~/public/Projekty/KS/OmicSelector/vignettes") # change it you to your working directory
#' # OmicSelector_OmicSelector(m = i, max_iterations = 1, stamp = "tutorial", debug = T) # we set debug to get more output
#' # }
#' # stopCluster(cl)
#'
#' @import snow doParallel
#'
#' @export
OmicSelector_OmicSelector = function(wd = getwd(), m = c(1:70),
                            max_iterations = 10, code_path = system.file("extdata", "", package = "OmicSelector"),
                            register_parallel = T, clx = NULL, stamp = as.numeric(Sys.time()),
                            prefer_no_features = 11, conda_path = "/home/konrad/anaconda3/bin/conda", debug = F,
                            timeout_sec = 172800, type = "auto") {

  oldwd = getwd()
  setwd(wd)
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
  suppressMessages(library(R.utils))
  suppressMessages(library(doParallel))


  if(!dir.exists("temp")) { dir.create("temp") }

  run_id = stamp
  formulas = list()
  times = list()

  zz <- file(paste0("temp/",stamp,paste0(m, collapse = "+"),"featureselection.log"), open = "wt")
  sink(zz)
  #sink(zz, type = "message")
  #pdf(paste0("temp/",stamp,paste0(m, collapse = "+"),"featureselection.pdf"))

  # not_needed as set to auto
  # if(file.exists("var_type.txt")) { type = readLines("var_type.txt", warn = F) } 

  wynik_finalny = withTimeout({
  dane = OmicSelector_load_datamix(); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]

  # n = 1
  n= 1
  if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Staring method n = ", n, ". Processing...")); OmicSelector_log(logfile = "temp/featureselection.log");
    start_time <- Sys.time()

    formulas[["all"]] = OmicSelector_create_formula(colnames(trainx))

    end_time <- Sys.time()
    saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS"))
    saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
    if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Standard DE...")
  
  wyniki = OmicSelector_differential_expression_ttest(trainx, train$Class, mode = type)
  istotne = filter(wyniki, `p-value` <= 0.05) %>% arrange(`p-value`)
  if(length(istotne$miR) == 0) { istotne = wyniki %>% arrange(`p-value`) } # what to do if none significant
  istotne_top = wyniki %>% arrange(`p-value`) %>% head(prefer_no_features)
  istotne_topBonf = wyniki %>% arrange(`p-value Bonferroni`) %>% head(prefer_no_features)
  istotne_topHolm = wyniki %>% arrange(`p-value Holm`) %>% head(prefer_no_features)
  istotne_topFC = wyniki %>% arrange(desc(abs(`log2FC`))) %>% head(prefer_no_features)

  train_sig = dplyr::select(train, as.character(istotne$miR), Class)
  trainx_sig = dplyr::select(train, as.character(istotne$miR))

  OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Preparing for SMOTE...")
  #train_sig_smoted = DMwR::SMOTE(Class ~ ., data = train_sig, perc.over = 10000,perc.under=100, k=10)
  train_sig_smoted = dplyr::select(train_smoted, as.character(istotne$miR), Class)
  train_sig_smoted$Class = factor(train_sig_smoted$Class, levels = c("Control","Case"))
  trainx_sig_smoted = dplyr::select(train_sig_smoted, starts_with("hsa"))

  wyniki_smoted = OmicSelector_differential_expression_ttest(trainx_smoted, train_smoted$Class, mode = type)
  istotne_smoted = filter(wyniki_smoted, `p-value BH` <= 0.05) %>% arrange(`p-value BH`)
  if(length(istotne_smoted$miR) == 0) { istotne = wyniki_smoted %>% arrange(`p-value BH`) } # what to do if none significant
  istotne_top_smoted = wyniki_smoted %>% arrange(`p-value BH`) %>% head(prefer_no_features)
  istotne_topBonf_smoted = wyniki_smoted %>% arrange(`p-value Bonferroni`) %>% head(prefer_no_features)
  istotne_topHolm_smoted = wyniki_smoted %>% arrange(`p-value Holm`) %>% head(prefer_no_features)
  istotne_topFC_smoted = wyniki_smoted %>% arrange(desc(abs(`log2FC`))) %>% head(prefer_no_features)

  # Caret prep
  suppressMessages(library(doParallel))
  registerDoSEQ()
  if (register_parallel) {
    OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Getting cluster ready for parallel computations...")
    if(is.null(clx)) {
      cl <- makePSOCKcluster(useXDR = TRUE, detectCores())
      registerDoParallel(cl)
      #on.exit(stopCluster(cl))
      }
    else { registerDoParallel(clx)
    #on.exit(stopCluster(clx))
      }
  }

  # 0. All and sig
  OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Cluster prepared. Moving to method matching...")

  # n = 2
  n = n + 1
  if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting.."));
    OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting SIG")
    start_time <- Sys.time()

    formulas[["sig"]] = OmicSelector_create_formula(as.character(istotne$miR))
    formulas[["sigtop"]] = OmicSelector_create_formula(as.character(istotne_top$miR))
    formulas[["sigtopBonf"]] = OmicSelector_create_formula(as.character(istotne_topBonf$miR))
    formulas[["sigtopHolm"]] = OmicSelector_create_formula(as.character(istotne_topHolm$miR))
    formulas[["topFC"]] = OmicSelector_create_formula(as.character(istotne_topFC$miR))
    formulas[["sigSMOTE"]] = OmicSelector_create_formula(as.character(istotne_smoted$miR))
    formulas[["sigtopSMOTE"]] = OmicSelector_create_formula(as.character(istotne_top_smoted$miR))
    formulas[["sigtopBonfSMOTE"]] = OmicSelector_create_formula(as.character(istotne_topBonf_smoted$miR))
    formulas[["sigtopHolmSMOTE"]] = OmicSelector_create_formula(as.character(istotne_topHolm_smoted$miR))
    formulas[["topFCSMOTE"]] = OmicSelector_create_formula(as.character(istotne_topFC_smoted$miR))


    end_time <- Sys.time()
    saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS"))
    saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
    if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 3
  # 1. Fold-change and sig. filter
  n = n + 1
  if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting.."));
    start_time <- Sys.time()

    OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting FC and SIG")
    fcsig = as.character(istotne$miR[abs(istotne$log2FC)>1])
    formulas[["fcsig"]] = OmicSelector_create_formula(fcsig)

    fcsig = as.character(istotne_smoted$miR[abs(istotne_smoted$log2FC)>1])
    formulas[["fcsigSMOTE"]] = OmicSelector_create_formula(fcsig)

    end_time <- Sys.time()
    saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS"))
    saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))

    if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 4
  # 2. CFS
  n = n + 1
  if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting.."));
    start_time <- Sys.time()
    OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting CFS")
    cfs = select.cfs(train)
    formulas[["cfs"]] = OmicSelector_create_formula(as.character(cfs$Biomarker))

    cfs = select.cfs(train_smoted)
    formulas[["cfsSMOTE"]] = OmicSelector_create_formula(as.character(cfs$Biomarker))

    cfs = select.cfs(train_sig)
    formulas[["cfs_sig"]] = OmicSelector_create_formula(as.character(cfs$Biomarker))

    cfs = select.cfs(train_sig_smoted)
    formulas[["cfsSMOTE_sig"]] = OmicSelector_create_formula(as.character(cfs$Biomarker))
    end_time <- Sys.time()
    saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS"))
    saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
    if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 5
  # 3. classifier.loop
  n = n + 1
  if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting.."));
    start_time <- Sys.time()
    OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting classifier loop")
    classloop = classifier.loop(train, feature.selection = "auc", method.cross="fold-crossval", classifiers=c("svm","lda","rf","nsc"), no.feat=prefer_no_features)
    f_classloop = rownames(classloop$no.selected)[classloop$no.selected[,1]>0]
    formulas[["classloop"]] = OmicSelector_create_formula(f_classloop)
    end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
    if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 6
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  classloop = classifier.loop(train_smoted, feature.selection = "auc", method.cross="fold-crossval", classifiers=c("svm","lda","rf","nsc"), no.feat=prefer_no_features)
  f_classloop = rownames(classloop$no.selected)[classloop$no.selected[,1]>0]
  formulas[["classloopSMOTE"]] = OmicSelector_create_formula(f_classloop)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 7
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  classloop = classifier.loop(train_sig, feature.selection = "auc", method.cross="fold-crossval", classifiers=c("svm","lda","rf","nsc"), no.feat=prefer_no_features)
  f_classloop = rownames(classloop$no.selected)[classloop$no.selected[,1]>0]
  formulas[["classloop_sig"]] = OmicSelector_create_formula(f_classloop)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 8
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  classloop = classifier.loop(train_sig_smoted, feature.selection = "auc", method.cross="fold-crossval", classifiers=c("svm","lda","rf","nsc"), no.feat=prefer_no_features)
  f_classloop = rownames(classloop$no.selected)[classloop$no.selected[,1]>0]
  formulas[["classloopSMOTE_sig"]] = OmicSelector_create_formula(f_classloop)

  end_time <- Sys.time()
  saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS"))
  saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 9
  # 4. select.forward.Corr
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting select.forward.Corr")
  fcfs = select.forward.Corr(train, disc.method="MDL")
  formulas[["fcfs"]] = OmicSelector_create_formula(fcfs)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 10
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  fcfs = select.forward.Corr(train_smoted, disc.method="MDL")
  formulas[["fcfsSMOTE"]] = OmicSelector_create_formula(fcfs)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 11
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  fcfs = select.forward.Corr(train_sig, disc.method="MDL")
  formulas[["fcfs_sig"]] = OmicSelector_create_formula(fcfs)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 12
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  fcfs = select.forward.Corr(train_sig_smoted, disc.method="MDL")
  formulas[["fcfsSMOTE_sig"]] = OmicSelector_create_formula(fcfs)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 13
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  # 5. select.forward.wrapper
  fwrap = select.forward.wrapper(train)
  formulas[["fwrap"]] = OmicSelector_create_formula(fwrap)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 14
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  fwrap = select.forward.wrapper(train_smoted)
  formulas[["fwrapSMOTE"]] = OmicSelector_create_formula(fwrap)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 15
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  fwrap = select.forward.wrapper(train_sig)
  formulas[["fwrap_sig"]] = OmicSelector_create_formula(fwrap)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 16
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  fwrap = select.forward.wrapper(train_sig_smoted)
  formulas[["fwrapSMOTE_sig"]] = OmicSelector_create_formula(fwrap)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }
  # 6, 7, 8.
  #select.process(dattable,method="InformationGain",disc.method="MDL",
  #               threshold=0.2,threshold.consis=0.05,attrs.nominal=numeric(),
  #               max.no.features=10)
  # n = 17
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting MDL")
  formulas[["AUC_MDL"]] = OmicSelector_create_formula(colnames(train)[select.process(train, method="auc", disc.method = "MDL", max.no.features = prefer_no_features)])
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 18
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  formulas[["SU_MDL"]] = OmicSelector_create_formula(colnames(train)[select.process(train, method="symmetrical.uncertainty", disc.method = "MDL", max.no.features = prefer_no_features)])
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 19
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  formulas[["CorrSF_MDL"]] = OmicSelector_create_formula(colnames(train)[select.process(train, method="CorrSF", disc.method = "MDL", max.no.features = prefer_no_features)])
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 20
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();

  formulas[["AUC_MDLSMOTE"]] = OmicSelector_create_formula(colnames(train_smoted)[select.process(train_smoted, method="auc", disc.method = "MDL", max.no.features = prefer_no_features)])
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 21
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();


  formulas[["SU_MDLSMOTE"]] = OmicSelector_create_formula(colnames(train_smoted)[select.process(train_smoted, method="symmetrical.uncertainty", disc.method = "MDL", max.no.features = prefer_no_features)])
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 22
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  formulas[["CorrSF_MDLSMOTE"]] = OmicSelector_create_formula(colnames(train_smoted)[select.process(train_smoted, method="CorrSF", disc.method = "MDL", max.no.features = prefer_no_features)])
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 23
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  formulas[["AUC_MDL_sig"]] = OmicSelector_create_formula(colnames(train_sig)[select.process(train_sig, method="auc", disc.method = "MDL", max.no.features = prefer_no_features)])
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 24
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  formulas[["SU_MDL_sig"]] = OmicSelector_create_formula(colnames(train_sig)[select.process(train_sig, method="symmetrical.uncertainty", disc.method = "MDL", max.no.features = prefer_no_features)])

  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 25
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  formulas[["CorrSF_MDL_sig"]] = OmicSelector_create_formula(colnames(train_sig)[select.process(train_sig, method="CorrSF", disc.method = "MDL", max.no.features = prefer_no_features)])
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 26
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  formulas[["AUC_MDLSMOTE_sig"]] = OmicSelector_create_formula(colnames(train_sig_smoted)[select.process(train_sig_smoted, method="auc", disc.method = "MDL", max.no.features = prefer_no_features)])
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 27
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  formulas[["SU_MDLSMOTE_sig"]] = OmicSelector_create_formula(colnames(train_sig_smoted)[select.process(train_sig_smoted, method="symmetrical.uncertainty", disc.method = "MDL", max.no.features = prefer_no_features)])
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 28
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  formulas[["CorrSF_MDLSMOTE_sig"]] = OmicSelector_create_formula(colnames(train_sig_smoted)[select.process(train_sig_smoted, method="CorrSF", disc.method = "MDL", max.no.features = prefer_no_features)])
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 29
  # 9. bounceR - genetic
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting bounceR..")
  suppressMessages(library(bounceR))
  mrmr <- bounceR::featureSelection(data = train,
                                    target = "Class",
                                    max_time = "15 mins",
                                    selection = selectionControl(n_rounds = NULL,
                                                                 n_mods = NULL,
                                                                 p = prefer_no_features,
                                                                 penalty = 0.5,
                                                                 reward = 0.2),
                                    bootstrap = "regular",
                                    early_stopping = "aic",
                                    boosting = boostingControl(mstop = 100, nu = 0.1),
                                    cores = parallel::detectCores()-1)
  formulas[["bounceR-full"]] = mrmr@opt_formula
  formulas[["bounceR-stability"]] = OmicSelector_create_formula(as.character(mrmr@stability[1:prefer_no_features,] %>% pull('feature')))
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 30
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  suppressMessages(library(bounceR))
  mrmr <- bounceR::featureSelection(data = train_smoted,
                           target = "Class",
                           max_time = "15 mins",
                           selection = selectionControl(n_rounds = NULL,
                                                        n_mods = NULL,
                                                        p = prefer_no_features,
                                                        penalty = 0.5,
                                                        reward = 0.2),
                           bootstrap = "regular",
                           early_stopping = "aic",
                           boosting = boostingControl(mstop = 100, nu = 0.1),
                           cores = parallel::detectCores()-1)
  formulas[["bounceR-full_SMOTE"]] = mrmr@opt_formula
  formulas[["bounceR-stability_SMOTE"]] = OmicSelector_create_formula(as.character(mrmr@stability[1:prefer_no_features,] %>% pull('feature')))
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 31
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  suppressMessages(library(bounceR))
  mrmr <- bounceR::featureSelection(data = train_sig,
                           target = "Class",
                           max_time = "15 mins",
                           selection = selectionControl(n_rounds = NULL,
                                                        n_mods = NULL,
                                                        p = prefer_no_features,
                                                        penalty = 0.5,
                                                        reward = 0.2),
                           bootstrap = "regular",
                           early_stopping = "aic",
                           boosting = boostingControl(mstop = 100, nu = 0.1),
                           cores = parallel::detectCores()-1)
  formulas[["bounceR-full_SIG"]] = mrmr@opt_formula
  formulas[["bounceR-stability_SIG"]] = OmicSelector_create_formula(as.character(mrmr@stability[1:prefer_no_features,] %>% pull('feature')))
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 32
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  suppressMessages(library(bounceR))
  mrmr <- bounceR::featureSelection(data = train_sig_smoted,
                           target = "Class",
                           max_time = "15 mins",
                           selection = selectionControl(n_rounds = NULL,
                                                        n_mods = NULL,
                                                        p = prefer_no_features,
                                                        penalty = 0.5,
                                                        reward = 0.2),
                           bootstrap = "regular",
                           early_stopping = "aic",
                           boosting = boostingControl(mstop = 100, nu = 0.1),
                           cores = parallel::detectCores()-1)
  formulas[["bounceR-full_SIGSMOTE"]] = mrmr@opt_formula
  formulas[["bounceR-stability_SIGSMOTE"]] = OmicSelector_create_formula(as.character(mrmr@stability[1:prefer_no_features,] %>% pull('feature')))
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 33
  # 10. RFE RandomForest
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting RFE RF")
  ctrl <- rfeControl(functions =rfFuncs,
                     method = "cv", number = 10,
                     saveDetails = TRUE,
                     allowParallel = TRUE,
                     returnResamp = "all",
                     verbose = T)

  rfProfile <- rfe(trainx, train$Class,
                   sizes = 3:11,
                   rfeControl = ctrl)
  plot(rfProfile, type=c("g", "o"))
  formulas[["RandomForestRFE"]] = OmicSelector_create_formula(predictors(rfProfile))
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 34
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  ctrl <- rfeControl(functions =rfFuncs,
                     method = "cv", number = 10,
                     saveDetails = TRUE,
                     allowParallel = TRUE,
                     returnResamp = "all",
                     verbose = T)
  rfProfile <- rfe(trainx_smoted, train_smoted$Class,
                   sizes = 3:11,
                   rfeControl = ctrl)
  plot(rfProfile, type=c("g", "o"))
  formulas[["RandomForestRFESMOTE"]] = OmicSelector_create_formula(predictors(rfProfile))
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 35
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  ctrl <- rfeControl(functions =rfFuncs,
                     method = "cv", number = 10,
                     saveDetails = TRUE,
                     allowParallel = TRUE,
                     returnResamp = "all",
                     verbose = T)
  rfProfile <- rfe(trainx_sig, train_sig$Class,
                   sizes = 3:11,
                   rfeControl = ctrl)
  plot(rfProfile, type=c("g", "o"))
  formulas[["RandomForestRFE_sig"]] = OmicSelector_create_formula(predictors(rfProfile))
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 36
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  ctrl <- rfeControl(functions =rfFuncs,
                     method = "cv", number = 10,
                     saveDetails = TRUE,
                     allowParallel = TRUE,
                     returnResamp = "all",
                     verbose = T)
  rfProfile <- rfe(trainx_sig_smoted, train_sig_smoted$Class,
                   sizes = 3:11,
                   rfeControl = ctrl)
  plot(rfProfile, type=c("g", "o"))
  formulas[["RandomForestRFESMOTE_sig"]] = OmicSelector_create_formula(predictors(rfProfile))
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }


  # n = 37
  # 11. Genetic
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting Genetic")
  ga_ctrl <- gafsControl(functions = rfGA, method = "repeatedcv", number=10, repeats=5, allowParallel=T)
  rf_ga <- gafs(x = trainx, y = train$Class,
                iters = max_iterations,
                gafsControl = ga_ctrl)
  plot(rf_ga) + theme_bw()
  print(rf_ga)
  formulas[["GeneticAlgorithmRF"]] = OmicSelector_create_formula(rf_ga$ga$final)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 38
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  ga_ctrl <- gafsControl(functions = rfGA, method = "repeatedcv", number=10, repeats=5, allowParallel=T)
  rf_ga <- gafs(x = trainx_smoted, y = train_smoted$Class,
                iters = max_iterations,
                gafsControl = ga_ctrl)
  plot(rf_ga) + theme_bw()
  print(rf_ga)
  formulas[["GeneticAlgorithmRFSMOTE"]] = OmicSelector_create_formula(rf_ga$ga$final)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 39
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  ga_ctrl <- gafsControl(functions = rfGA, method = "repeatedcv", number=10, repeats=5, allowParallel=T)
  rf_ga <- gafs(x = trainx_sig, y = train_sig$Class,
                iters = max_iterations,
                gafsControl = ga_ctrl)
  plot(rf_ga) + theme_bw()
  print(rf_ga)
  formulas[["GeneticAlgorithmRF_sig"]] = OmicSelector_create_formula(rf_ga$ga$final)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 40
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  ga_ctrl <- gafsControl(functions = rfGA, method = "repeatedcv", number=10, repeats=5, allowParallel=T)
  rf_ga <- gafs(x = trainx_sig_smoted, y = train_sig_smoted$Class,
                iters = max_iterations,
                gafsControl = ga_ctrl)
  plot(rf_ga) + theme_bw()
  print(rf_ga)
  formulas[["GeneticAlgorithmRFSMOTE_sig"]] = OmicSelector_create_formula(rf_ga$ga$final)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 41
  # 12. SimulatedAnealing
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting SimulatedAnealing")
  sa_ctrl <- safsControl(functions = rfSA,
                         method = "repeatedcv",
                         number=5, repeats=10, allowParallel=T,
                         improve = 50)

  rf_sa <- safs(x = trainx, y = train$Class,
                iters = max_iterations,
                safsControl = sa_ctrl)
  print(rf_sa)
  plot(rf_sa) + theme_bw()
  formulas[["SimulatedAnnealingRF"]] = OmicSelector_create_formula(rf_sa$sa$final)



  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 42
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
    OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting SimulatedAnealing")
  sa_ctrl <- safsControl(functions = rfSA,
                         method = "repeatedcv",
                         number=5, repeats=10, allowParallel=T,
                         improve = 50)
  
  rf_sa <- safs(x = trainx_smoted, y = train_smoted$Class,
                iters = max_iterations,
                safsControl = sa_ctrl)
  print(rf_sa)
  plot(rf_sa) + theme_bw()
  formulas[["SimulatedAnnealingRFSMOTE"]] = OmicSelector_create_formula(rf_sa$sa$final)

  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 43
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
    OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting SimulatedAnealing")
  sa_ctrl <- safsControl(functions = rfSA,
                         method = "repeatedcv",
                         number=5, repeats=10, allowParallel=T,
                         improve = 50)
  
  rf_sa <- safs(x = trainx_sig, y = train_sig$Class,
                iters = max_iterations,
                safsControl = sa_ctrl)
  print(rf_sa)
  plot(rf_sa) + theme_bw()
  formulas[["SimulatedAnnealingRF_sig"]] = OmicSelector_create_formula(rf_sa$sa$final)

  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 44
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
    OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting SimulatedAnealing")
  sa_ctrl <- safsControl(functions = rfSA,
                         method = "repeatedcv",
                         number=5, repeats=10, allowParallel=T,
                         improve = 50)
  
  rf_sa <- safs(x = trainx_sig_smoted, y = train_sig_smoted$Class,
                iters = max_iterations,
                safsControl = sa_ctrl)
  print(rf_sa)
  plot(rf_sa) + theme_bw()
  formulas[["SimulatedAnnealingRFSMOTE_sig"]] = OmicSelector_create_formula(rf_sa$sa$final)

  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }


  # n = 45
  # 14. Boruta (https://www.jstatsoft.org/article/view/v036i11)
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  suppressMessages(library(Boruta))
  OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting Boruta")
  bor = Boruta(trainx, train$Class)
  formulas[["Boruta"]] = OmicSelector_create_formula(names(bor$finalDecision)[as.character(bor$finalDecision) == "Confirmed"])

  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 46
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  bor = Boruta(trainx_smoted, train_smoted$Class)
  formulas[["BorutaSMOTE"]] = OmicSelector_create_formula(names(bor$finalDecision)[as.character(bor$finalDecision) == "Confirmed"])

  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 47
  # 15. spFSR - feature selection and ranking by simultaneous perturbation stochastic approximation
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting spFSR")
  suppressMessages(library(spFSR))



  knnWrapper    <- makeLearner("classif.knn", k = 5)
  classifTask   <- makeClassifTask(data = train, target = "Class")
  perf.measure  <- acc
  spsaMod <- spFeatureSelection(
    task = classifTask,
    wrapper = knnWrapper,
    measure = perf.measure ,
    num.features.selected = prefer_no_features,
    iters.max = max_iterations,
    num.cores = detectCores())
  formulas[["spFSR"]] = OmicSelector_create_formula(spsaMod$features)


  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 48
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting spFSR")
  suppressMessages(library(spFSR))


  
  classifTask   <- makeClassifTask(data = train_smoted, target = "Class")
  perf.measure  <- acc
  spsaMod <- spFeatureSelection(
    task = classifTask,
    wrapper = knnWrapper,
    measure = perf.measure ,
    num.features.selected = prefer_no_features,
    iters.max = max_iterations,
    num.cores = detectCores())
  formulas[["spFSRSMOTE"]] = OmicSelector_create_formula(spsaMod$features)

  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 49
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  # varSelRF
  OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting varSelRF")
  suppressMessages(library(varSelRF))
  var.sel <- varSelRF(trainx, train$Class, ntree = 500, ntreeIterat = max_iterations, vars.drop.frac = 0.05, whole.range = T, keep.forest = T)
  formulas[["varSelRF"]] = OmicSelector_create_formula(var.sel$selected.vars)

  var.sel <- varSelRF(trainx_smoted, train_smoted$Class, ntree = 500, ntreeIterat = max_iterations, vars.drop.frac = 0.05, whole.range = T, keep.forest = T)
  formulas[["varSelRFSMOTE"]] = OmicSelector_create_formula(var.sel$selected.vars)

  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 50
  n = n + 1; if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  # 13. WxNet (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6642261/)
  OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting WxNet")
  suppressMessages(library(plyr))
  suppressMessages(library(dplyr))
  suppressMessages(library(reticulate))
  suppressMessages(library(tidyverse))
  suppressMessages(library(data.table))
  suppressMessages(library(DMwR))



  # Set the path to the Python executale file
  #use_python("/anaconda3/bin/python", required = T)


  # conda_list()
  try({ conda_create("wxnet", c("tensorflow-gpu","keras")) })
  use_condaenv("wxnet", required = T)
  py_config()

  dane = OmicSelector_load_datamix(replace_smote = F); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]

  #train_org = train
  #train = SMOTE(Class ~ ., data = train_org, perc.over = 10000,perc.under=100)

  # Przygotowanie
  trainx = dplyr::select(train, starts_with("hsa"))
  testx = dplyr::select(test, starts_with("hsa"))
  trainx = t(trainx)
  testx = t(testx)
  ids = paste0("train",1:ncol(trainx))
  ids2 = paste0("test",1:ncol(testx))
  colnames(trainx) = gsub("-", "", ids)
  colnames(testx) = gsub("-", "", ids2)

  traindata = cbind(data.frame(fnames = rownames(trainx), trainx))
  fwrite(traindata, paste0(code_path, "wx/DearWXpub/src/train-data.csv"), row.names = F, quote = F)
  testdata = cbind(data.frame(fnames = rownames(testx), testx))
  fwrite(testdata, paste0(code_path, "wx/DearWXpub/src/test-data.csv"), row.names = F, quote = F)

  trainanno = data.frame(id = colnames(trainx), label = train$Class)
  fwrite(trainanno, paste0(code_path, "wx/DearWXpub/src/train-anno.csv"), row.names = F, quote = F)
  testanno = data.frame(id = colnames(testx), label = test$Class)
  fwrite(testanno, paste0(code_path, "wx/DearWXpub/src/test-anno.csv"), row.names = F, quote = F)


  # Wywolanie WX

  out <- tryCatch(
    {
      setwd(paste0(code_path,"wx/DearWXpub/src/"))
      # system(paste0(conda_path," activate base"))
      py_run_file("wx_konsta.py", local = T)

      np <- import("numpy")
      w = np$load("wyniki.npy",allow_pickle=T)
      formulas[["Wx"]] = OmicSelector_create_formula(w)
    },
    error=function(cond) {
      message("ERROR:")
      message(cond)
      # Choose a return value in case of error
    },
    warning=function(cond) {
      message("WARNING:")
      message(cond)
    },
    finally={
      setwd(wd)
    }
  )





  # Wx with SMOTE
  dane = OmicSelector_load_datamix(replace_smote = T); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]

  # Przygotowanie
  trainx = dplyr::select(train, starts_with("hsa"))
  testx = dplyr::select(test, starts_with("hsa"))
  trainx = t(trainx)
  testx = t(testx)
  ids = paste0("train",1:ncol(trainx))
  ids2 = paste0("test",1:ncol(testx))
  colnames(trainx) = gsub("-", "", ids)
  colnames(testx) = gsub("-", "", ids2)

  traindata = cbind(data.frame(fnames = rownames(trainx), trainx))
  fwrite(traindata, paste0(code_path, "wx/DearWXpub/src/train-data.csv"), row.names = F, quote = F)
  testdata = cbind(data.frame(fnames = rownames(testx), testx))
  fwrite(testdata, paste0(code_path, "wx/DearWXpub/src/test-data.csv"), row.names = F, quote = F)

  trainanno = data.frame(id = colnames(trainx), label = train$Class)
  fwrite(trainanno, paste0(code_path, "wx/DearWXpub/src/train-anno.csv"), row.names = F, quote = F)
  testanno = data.frame(id = colnames(testx), label = test$Class)
  fwrite(testanno, paste0(code_path, "wx/DearWXpub/src/test-anno.csv"), row.names = F, quote = F)


  # Wywolanie WX
  # setwd(paste0(code_path,"wx/DearWXpub/src/"))
  # py_run_file("wx_konsta.py", local = T)
  #
  # np <- import("numpy")
  # w = np$load("wyniki.npy",allow_pickle=T)
  # print(w)
  # setwd(wd)
  # formulas[["WxSMOTE"]] = OmicSelector_create_formula(w)
  out <- tryCatch(
    {
      setwd(paste0(code_path,"wx/DearWXpub/src/"))
      py_run_file("wx_konsta.py", local = T)

      np <- import("numpy")
      w = np$load("wyniki.npy",allow_pickle=T)
      formulas[["WxSMOTE"]] = OmicSelector_create_formula(w)
    },
    error=function(cond) {
      message("ERROR:")
      message(cond)
      # Choose a return value in case of error
    },
    warning=function(cond) {
      message("WARNING:")
      message(cond)
    },
    finally={
      setwd(wd)
    }
  )
  # Salvage results if SMOTE fails
  saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))

  dane = OmicSelector_load_datamix(replace_smote = F); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]

  # Przygotowanie
  trainx = dplyr::select(train, starts_with("hsa"))
  testx = dplyr::select(test, starts_with("hsa"))
  trainx = t(trainx)
  testx = t(testx)
  ids = paste0("train",1:ncol(trainx))
  ids2 = paste0("test",1:ncol(testx))
  colnames(trainx) = gsub("-", "", ids)
  colnames(testx) = gsub("-", "", ids2)

  traindata = cbind(data.frame(fnames = rownames(trainx), trainx))
  fwrite(traindata, paste0(code_path, "wx/DearWXpub/src/train-data.csv"), row.names = F, quote = F)
  testdata = cbind(data.frame(fnames = rownames(testx), testx))
  fwrite(testdata, paste0(code_path, "wx/DearWXpub/src/test-data.csv"), row.names = F, quote = F)

  trainanno = data.frame(id = colnames(trainx), label = train$Class)
  fwrite(trainanno, paste0(code_path, "wx/DearWXpub/src/train-anno.csv"), row.names = F, quote = F)
  testanno = data.frame(id = colnames(testx), label = test$Class)
  fwrite(testanno, paste0(code_path, "wx/DearWXpub/src/test-anno.csv"), row.names = F, quote = F)


  # Wywolanie WX
  # setwd(paste0(code_path,"wx/DearWXpub/src/"))
  # py_run_file("wx_konsta_z.py", local = T)
  #
  # np <- import("numpy")
  # w = np$load("wyniki.npy",allow_pickle=T)
  # print(w)
  # setwd(wd)
  # formulas[["Wx_Zscore"]] = OmicSelector_create_formula(w)
  out <- tryCatch(
    {
      setwd(paste0(code_path,"wx/DearWXpub/src/"))
      py_run_file("wx_konsta_z.py", local = T)

      np <- import("numpy")
      w = np$load("wyniki.npy",allow_pickle=T)
      formulas[["Wx_Zscore"]] = OmicSelector_create_formula(w)
    },
    error=function(cond) {
      message("ERROR:")
      message(cond)
      # Choose a return value in case of error
    },
    warning=function(cond) {
      message("WARNING:")
      message(cond)
    },
    finally={
      setwd(wd)
    }
  )

  dane = OmicSelector_load_datamix(replace_smote = T); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]

  # Przygotowanie
  trainx = dplyr::select(train, starts_with("hsa"))
  testx = dplyr::select(test, starts_with("hsa"))
  trainx = t(trainx)
  testx = t(testx)
  ids = paste0("train",1:ncol(trainx))
  ids2 = paste0("test",1:ncol(testx))
  colnames(trainx) = gsub("-", "", ids)
  colnames(testx) = gsub("-", "", ids2)

  traindata = cbind(data.frame(fnames = rownames(trainx), trainx))
  fwrite(traindata, paste0(code_path, "wx/DearWXpub/src/train-data.csv"), row.names = F, quote = F)
  testdata = cbind(data.frame(fnames = rownames(testx), testx))
  fwrite(testdata, paste0(code_path, "wx/DearWXpub/src/test-data.csv"), row.names = F, quote = F)

  trainanno = data.frame(id = colnames(trainx), label = train$Class)
  fwrite(trainanno, paste0(code_path, "wx/DearWXpub/src/train-anno.csv"), row.names = F, quote = F)
  testanno = data.frame(id = colnames(testx), label = test$Class)
  fwrite(testanno, paste0(code_path, "wx/DearWXpub/src/test-anno.csv"), row.names = F, quote = F)


  # Wywolanie WX
  # setwd(paste0(code_path,"wx/DearWXpub/src/"))
  # py_run_file("wx_konsta_z.py", local = T)
  #
  # np <- import("numpy")
  # w = np$load("wyniki.npy",allow_pickle=T)
  # print(w)
  # setwd(wd)
  # formulas[["Wx_ZscoreSMOTE"]] = OmicSelector_create_formula(w)
  out <- tryCatch(
    {
      setwd(paste0(code_path,"wx/DearWXpub/src/"))
      py_run_file("wx_konsta_z.py", local = T)

      np <- import("numpy")
      w = np$load("wyniki.npy",allow_pickle=T)
      formulas[["Wx_ZscoreSMOTE"]] = OmicSelector_create_formula(w)
    },
    error=function(cond) {
      message("ERROR:")
      message(cond)
      # Choose a return value in case of error
    },
    warning=function(cond) {
      message("WARNING:")
      message(cond)
    },
    finally={
      setwd(wd)
    }
  )


  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 51
  # My.stepwise
  n = n + 1
  if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting.."));
    OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting My.stepwise")
    start_time <- Sys.time()

    dane = OmicSelector_load_datamix(replace_smote = F); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]

    suppressMessages(library(My.stepwise))
    temp = capture.output(My.stepwise.glm(Y = "Class", colnames(trainx), data = train, sle = 0.05, sls = 0.05, myfamily = "binomial"))
    # temp2 = temp[length(temp)-1]
    # temp3 = temp[length(temp)-3]
    # temp4 = temp[length(temp)-5]
    temp5 = paste(temp[(length(temp)-12):length(temp)], collapse = " ")
    wybrane = FALSE
    for(i in 1:length(colnames(trainx)))
    {
      temp6 = colnames(trainx)[i]
      wybrane[i] = grepl(temp6, temp5)
    }
    formulas[["Mystepwise_glm_binomial"]] = OmicSelector_create_formula(colnames(trainx)[wybrane])

    wyniki = OmicSelector_differential_expression_ttest(trainx, train$Class, mode = type)
    istotne = filter(wyniki, `p-value BH` <= 0.05) %>% arrange(`p-value BH`)

    temp = capture.output(My.stepwise.glm(Y = "Class", as.character(istotne$miR), data = train, sle = 0.05, sls = 0.05, myfamily = "binomial"))
    # temp2 = temp[length(temp)-1]
    # temp3 = temp[length(temp)-3]
    # temp4 = temp[length(temp)-5]
    temp5 = paste(temp[(length(temp)-12):length(temp)], collapse = " ")
    wybrane = FALSE
    for(i in 1:length(colnames(trainx)))
    {
      temp6 = colnames(trainx)[i]
      wybrane[i] = grepl(temp6, temp5)
    }
    formulas[["Mystepwise_sig_glm_binomial"]] = OmicSelector_create_formula(colnames(trainx)[wybrane])

    end_time <- Sys.time()
    saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS"))
    saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
    if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 52
  n = n + 1
  if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting.."));
    OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting My.stepwise SMOTE")
    start_time <- Sys.time()

    dane = OmicSelector_load_datamix(replace_smote = F); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]

    suppressMessages(library(My.stepwise))
    temp = capture.output(My.stepwise.glm(Y = "Class", colnames(trainx_smoted), data = train_smoted, sle = 0.05, sls = 0.05, myfamily = "binomial"))
    # temp2 = temp[length(temp)-1]
    # temp3 = temp[length(temp)-3]
    # temp4 = temp[length(temp)-5]
    temp5 = paste(temp[(length(temp)-12):length(temp)], collapse = " ")
    wybrane = FALSE
    for(i in 1:length(colnames(trainx_smoted)))
    {
      temp6 = colnames(trainx_smoted)[i]
      wybrane[i] = grepl(temp6, temp5)
    }
    formulas[["Mystepwise_glm_binomialSMOTE"]] = OmicSelector_create_formula(colnames(trainx_smoted)[wybrane])

    wyniki = OmicSelector_differential_expression_ttest(trainx, train$Class, mode = type)
    istotne = filter(wyniki, `p-value BH` <= 0.05) %>% arrange(`p-value BH`)

    temp = capture.output(My.stepwise.glm(Y = "Class", as.character(istotne$miR), data = train_smoted, sle = 0.05, sls = 0.05, myfamily = "binomial"))
    # temp2 = temp[length(temp)-1]
    # temp3 = temp[length(temp)-3]
    # temp4 = temp[length(temp)-5]
    temp5 = paste(temp[(length(temp)-12):length(temp)], collapse = " ")
    wybrane = FALSE
    for(i in 1:length(colnames(trainx_smoted)))
    {
      temp6 = colnames(trainx_smoted)[i]
      wybrane[i] = grepl(temp6, temp5)
    }
    formulas[["Mystepwise_sig_glm_binomialSMOTE"]] = OmicSelector_create_formula(colnames(trainx_smoted)[wybrane])

    end_time <- Sys.time()
    saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS"))
    saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
    if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 53
  n = n + 1
  if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting.."));
    OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting stepAIC")
    start_time <- Sys.time()

    dane = OmicSelector_load_datamix(replace_smote = F); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]

    temp = glm(Class ~ ., data = train, family = "binomial")
    temp2 = stepAIC(temp)

    formulas[["stepAIC"]] = temp2$formula

    wyniki = OmicSelector_differential_expression_ttest(trainx, train$Class, mode = type)
    istotne = filter(wyniki, `p-value BH` <= 0.05) %>% arrange(`p-value BH`)

    train.sig = dplyr::select(train, as.character(istotne$miR), Class)

    temp = glm(Class ~ ., data = train.sig, family = "binomial")
    temp2 = stepAIC(temp)

    formulas[["stepAICsig"]] = temp2$formula

    end_time <- Sys.time()
    saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS"))
    saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
    if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 54
  n = n + 1
  if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting.."));
    OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting stepAIC SMOTE")
    start_time <- Sys.time()

    dane = OmicSelector_load_datamix(replace_smote = F); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]

    temp = glm(Class ~ ., data = train_smoted, family = "binomial")
    temp2 = stepAIC(temp)

    formulas[["stepAIC_SMOTE"]] = temp2$formula

    wyniki = OmicSelector_differential_expression_ttest(trainx, train$Class, mode = type)
    istotne = filter(wyniki, `p-value BH` <= 0.05) %>% arrange(`p-value BH`)

    train.sig = dplyr::select(train_smoted, as.character(istotne$miR), Class)

    temp = glm(Class ~ ., data = train.sig, family = "binomial")
    temp2 = stepAIC(temp)

    formulas[["stepAICsig_SMOTE"]] = temp2$formula

    end_time <- Sys.time()
    saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS"))
    saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
    if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 55
  n = n + 1
  if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting.."));
    OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting MK method (iterated RFE)")
    start_time <- Sys.time()

    dane = OmicSelector_load_datamix(replace_smote = F); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]

    selectedMirsCV <- OmicSelector_iteratedRFE(trainSet = train, useCV = T, classLab = 'Class', checkNFeatures = prefer_no_features)$topFeaturesPerN[[prefer_no_features]]
    selectedMirsTest <- mk.iteratedRFE(trainSet = train, testSet = test, classLab = 'Class', checkNFeatures = prefer_no_features)$topFeaturesPerN[[prefer_no_features]]

    formulas[["iteratedRFECV"]] = OmicSelector_create_formula(selectedMirsCV$topFeaturesPerN[[prefer_no_features]])
    formulas[["iteratedRFETest"]] = OmicSelector_create_formula(selectedMirsTest$topFeaturesPerN[[prefer_no_features]])

    end_time <- Sys.time()
    saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS"))
    saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
    if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 56
  n = n + 1
  if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting.."));
    OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting MK method (iterated RFE) with SMOTE")
    start_time <- Sys.time()

    dane = OmicSelector_load_datamix(replace_smote = F); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]

    

    formulas[["iteratedRFECV_SMOTE"]] = OmicSelector_create_formula(selectedMirsCV$topFeaturesPerN[[prefer_no_features]])
    formulas[["iteratedRFETest_SMOTE"]] = OmicSelector_create_formula(selectedMirsTest$topFeaturesPerN[[prefer_no_features]])

    end_time <- Sys.time()
    saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS"))
    saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
    if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }


  # n = 57
  # LASSO model with alpha = 1 (Least Absolute Shrinkage and Selection Operator - a type of regularization method that penalizes with L1-norm.). The function cv.glmnet() is used to search for a regularization parameter, namely Lambda, that controls the penalty strength
  n = n + 1
  if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting.."));
    OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting LASSO with and without SMOTE")
    start_time <- Sys.time()

    dane = OmicSelector_load_datamix(replace_smote = F); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]

    library("glmnet")
    lasso_fit <- cv.glmnet(as.matrix(trainx), train$Class, family = "binomial", alpha = 1)
    plot(lasso_fit)
    coef <- predict(lasso_fit, s = "lambda.min", type = "nonzero")
    result <- data.frame(GENE = names(as.matrix(coef(lasso_fit, s = "lambda.min"))
                                [as.matrix(coef(lasso_fit, s = "lambda.min"))[,1]!=0, 1])[-1], 
                   SCORE = as.numeric(as.matrix(coef(lasso_fit, s = "lambda.min"))
                                      [as.matrix(coef(lasso_fit, 
                                                      s = "lambda.min"))[,1]!=0, 1])[-1])
    result <- result[order(-abs(result$SCORE)),]
    formulas[["LASSO"]] = OmicSelector_create_formula(as.character(result$GENE))
    
    lasso_fit <- cv.glmnet(as.matrix(trainx_smoted), train_smoted$Class, family = "binomial", alpha = 1)
    plot(lasso_fit)
    coef <- predict(lasso_fit, s = "lambda.min", type = "nonzero")
    result <- data.frame(GENE = names(as.matrix(coef(lasso_fit, s = "lambda.min"))
                                [as.matrix(coef(lasso_fit, s = "lambda.min"))[,1]!=0, 1])[-1], 
                   SCORE = as.numeric(as.matrix(coef(lasso_fit, s = "lambda.min"))
                                      [as.matrix(coef(lasso_fit, 
                                                      s = "lambda.min"))[,1]!=0, 1])[-1])
    result <- result[order(-abs(result$SCORE)),]
    formulas[["LASSO_SMOTE"]] = OmicSelector_create_formula(as.character(result$GENE))

    end_time <- Sys.time()
    saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS"))
    saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
    if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # n = 58
  # Elastic net with tuning the value of Alpha through a line search with the parallelism
  n = n + 1
  if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting.."));
    OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting ElasticNet with and without SMOTE")
    start_time <- Sys.time()

    dane = OmicSelector_load_datamix(replace_smote = F); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]

    library(foreach)
    library(glmnet)

    a <- seq(0.1, 0.9, 0.05)
    search <- foreach(i = a, .combine = rbind) %dopar% {
        library(OmicSelector)
        dane = OmicSelector_load_datamix(replace_smote = F); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]
        library(glmnet)
        cv <- cv.glmnet(as.matrix(trainx), train$Class, family = "binomial", nfold = 10, parallel = TRUE, alpha = i)
        data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
    }
    cv3 <- search[search$cvm == min(search$cvm), ]
    md3 <- glmnet(as.matrix(trainx), train$Class, family = "binomial", lambda = cv3$lambda.1se, alpha = cv3$alpha)
    result <- data.frame(GENE = names(as.matrix(coef(md3, s = "lambda.min"))
                                      [as.matrix(coef(md3, s = "lambda.min"))[,1]!=0, 1])[-1], 
                          SCORE = as.numeric(as.matrix(coef(md3, s = "lambda.min"))
                                            [as.matrix(coef(md3, 
                                                            s = "lambda.min"))[,1]!=0, 1])[-1])
    result <- result[order(-abs(result$SCORE)),]
    formulas[["ElasticNet"]] = OmicSelector_create_formula(as.character(result$GENE))


    # SMOTE:
    a <- seq(0.1, 0.9, 0.05)
    search <- foreach(i = a, .combine = rbind) %dopar% {
        library(OmicSelector)
        dane = OmicSelector_load_datamix(replace_smote = F); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]
        library(glmnet)
        cv <- cv.glmnet(as.matrix(trainx_smoted), train_smoted$Class, family = "binomial", nfold = 10, parallel = TRUE, alpha = i)
        data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
    }
    cv3 <- search[search$cvm == min(search$cvm), ]
    md3 <- glmnet(as.matrix(trainx_smoted), train_smoted$Class, family = "binomial", lambda = cv3$lambda.1se, alpha = cv3$alpha)
    result <- data.frame(GENE = names(as.matrix(coef(md3, s = "lambda.min"))
                                      [as.matrix(coef(md3, s = "lambda.min"))[,1]!=0, 1])[-1], 
                          SCORE = as.numeric(as.matrix(coef(md3, s = "lambda.min"))
                                            [as.matrix(coef(md3, 
                                                            s = "lambda.min"))[,1]!=0, 1])[-1])
    result <- result[order(-abs(result$SCORE)),]
    formulas[["ElasticNet_SMOTE"]] = OmicSelector_create_formula(as.character(result$GENE))
    
    

    end_time <- Sys.time()
    saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS"))
    saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
    if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

  # # n = 59
  # # Ridge regression sorted
  # n = n + 1
  # if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting.."));
  #   OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting Ridge")
  #   start_time <- Sys.time()

  #   dane = OmicSelector_load_datamix(replace_smote = F); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]

  #   library(glmnet)
  #   fit <- cv.glmnet(as.matrix(trainx), train$Class, family = "binomial", alpha = 0)
  #   formulas[["Ridge"]] = OmicSelector_create_formula(as.character(names(sort(abs(coef(fit)[-1,]), decreasing = T)[1:prefer_no_features])))

  #   fit <- cv.glmnet(as.matrix(trainx_smoted), train_smoted$Class, family = "binomial", alpha = 0)
  #   formulas[["Ridge_SMOTE"]] = OmicSelector_create_formula(as.character(names(sort(abs(coef(fit)[-1,]), decreasing = T)[1:prefer_no_features])))

  #   end_time <- Sys.time()
  #   saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS"))
  #   saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  #   if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  # }

  # n = 59
  # Stepwise LDA
  n = n + 1
  if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting.."));
    OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting stepLDA")
    start_time <- Sys.time()

    dane = OmicSelector_load_datamix(replace_smote = F); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]

    tempdb = cbind(`Class` = train$Class, trainx)
    library(caret)
    slda <- train(Class ~ ., data = tempdb,
              method = "stepLDA",
              trControl = trainControl(method = "loocv"))
    formulas[["stepLDA"]] = OmicSelector_create_formula(predictors(slda))

tempdb = cbind(`Class` = train_smoted$Class, trainx_smoted)
    library(caret)
    slda <- train(Class ~ ., data = tempdb,
              method = "stepLDA",
              trControl = trainControl(method = "loocv"))
    formulas[["stepLDA_SMOTE"]] = OmicSelector_create_formula(predictors(slda))

    end_time <- Sys.time()
    saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS"))
    saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
    if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }


  # n = 60
  # feseR
  n = n + 1
  if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting.."));
    OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting feseR")
    start_time <- Sys.time()

    dane = OmicSelector_load_datamix(replace_smote = F); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]

    library(feseR)
    features <- as.matrix(scale(trainx, center=TRUE, scale=TRUE))
    rownames(features) = 1:nrow(features)
    cls <- as.numeric(train$Class)-1
    try({ output <- filter.corr(features = features, class = cls, mincorr = 0.2); formulas[["feseR_filter.corr"]] = OmicSelector_create_formula(colnames(output)); }, silent = T)
    try({ output <- filter.gain.inf(features = features, class = cls, zero.gain.out = TRUE); formulas[["feseR_gain.inf"]] = OmicSelector_create_formula(colnames(output)); }, silent = T)
    try({ output <- filter.matrix.corr(features = features, maxcorr = 0.75); formulas[["feseR_matrix.corr"]] = OmicSelector_create_formula(colnames(output)); }, silent = T)
    try({ output <- combineFS(features = features, class = cls,
                        univariate = 'corr', mincorr = 0.2,
                        multivariate = 'mcorr', maxcorr = 0.75,
                        wrapper = 'rfe.rf', number.cv = 10, 
                        group.sizes = 1:prefer_no_features, 
                        extfolds = 20); formulas[["feseR_combineFS_RF"]] = OmicSelector_create_formula(output$opt.variables); }, silent = T)

    features <- as.matrix(scale(trainx_smoted, center=TRUE, scale=TRUE))
    rownames(features) = 1:nrow(features)
    cls <- as.numeric(train_smoted$Class)-1
    try({ output <- filter.corr(features = features, class = cls, mincorr = 0.2); formulas[["feseR_filter.corr_SMOTE"]] = OmicSelector_create_formula(colnames(output)); }, silent = T)
    try({ output <- filter.gain.inf(features = features, class = cls, zero.gain.out = TRUE); formulas[["feseR_gain.inf_SMOTE"]] = OmicSelector_create_formula(colnames(output)); }, silent = T)
    try({ output <- filter.matrix.corr(features = features, maxcorr = 0.75); formulas[["feseR_matrix.corr_SMOTE"]] = OmicSelector_create_formula(colnames(output)); }, silent = T)
    try({ output <- combineFS(features = features, class = cls,
                        univariate = 'corr', mincorr = 0.2,
                        multivariate = 'mcorr', maxcorr = 0.75,
                        wrapper = 'rfe.rf', number.cv = 10, 
                        group.sizes = 1:prefer_no_features, 
                        extfolds = 20); formulas[["feseR_combineFS_RF_SMOTE"]] = OmicSelector_create_formula(output$opt.variables); }, silent = T)


    end_time <- Sys.time()
    saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS"))
    saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
    if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }

   # n = 61
  # Asakura2020
  n = n + 1
  if (n %in% m) { OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Matched method ", n, " with those requested.. Starting.."));
    OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = "Starting feseR")
    start_time <- Sys.time()

    dane = OmicSelector_load_datamix(replace_smote = F); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]

    source(system.file("extdata", "Asakura2020/init.R", package = "readr"))

    # TO DO

    fisher_pivot_num <- 20
    fisher_turn_num <- 3
    cancer_type <- "SC"
    output_dir <- "output"
    include_miRNA <- c()
    include_miR_logical <- "and"
    exclude_miRNA <- c()

    model <- train
    test <- test

  if (!identical(colnames(model), colnames(test))) {
    stop("MODEL/TEST column name mismatched")
  }

  for (i in seq_len(ncol(model))) {
    miR2 <- unlist(strsplit(colnames(model)[i], ","))
    if (length(miR2) > 1) {
      colnames(model)[i] <- miR2[1]
      colnames(test)[i] <- miR2[1]
    }
  }

  for (miR in exclude_miRNA) {
    col <- charmatch(miR, colnames(model), 0)
    if (col > 0) {
      model <- model[-col]
      test <- test[-col]
    }
  }

  fisher_miR <- matrix(0, nrow = fisher_pivot_num, ncol = fisher_turn_num + 1)

  print_log(paste0("Fisher process turn : ", 1))

  cols <- c()
  cols_hash <- hash()
  miRs <- next_miR(model, fisher_pivot_num, cols, cols_hash)
  sortlist <- order(miRs$V2, decreasing = T)
  miRs <- miRs[sortlist, ]

  for (i in seq_len(nrow(fisher_miR))) {
    fisher_miR[i, 1] <- miRs$V1[i]
    fisher_miR[i, ncol(fisher_miR)] <- miRs$V2[i]
  }
  output_turn_fisher(model, test, 1, fisher_miR, process_name)

  for (i in seq(2, fisher_turn_num)) {
    print_log(paste0("Fisher process turn : ", i))
    k <- 0
    wk_fisher_miR <- matrix(0,
      nrow = fisher_pivot_num * ncol(model),
      ncol = ncol(fisher_miR)
    )
    for (j in seq_len(fisher_pivot_num)) {
      if (is.na(fisher_miR[j, 1])) break
      cols <- fisher_miR[j, 1:(i - 1)]
      miRs <- next_miR(model, fisher_pivot_num, cols, cols_hash)
      d <- k * ncol(model)
      for (l in seq_len(nrow(miRs))) {
        wk_fisher_miR[d + l, 1:(i - 1)] <- cols
        wk_fisher_miR[d + l, i] <- miRs$V1[l]
        wk_fisher_miR[d + l, ncol(fisher_miR)] <- miRs$V2[l]
      }
      k <- k + 1
    }
    sortlist <- order(wk_fisher_miR[, ncol(fisher_miR)], decreasing = T)
    wk_fisher_miR <- wk_fisher_miR[sortlist, ]
    wk_fisher_miR <- wk_fisher_miR[wk_fisher_miR[, ncol(fisher_miR)] > 0, ]
    for (j in seq_len(min(nrow(fisher_miR), nrow(wk_fisher_miR)))) {
      fisher_miR[j, ] <- wk_fisher_miR[j, ]
    }
    output_turn_fisher(model, test, i, fisher_miR, process_name)
    clear(cols_hash)
  }


    


    end_time <- Sys.time()
    saveRDS(end_time - start_time, paste0("temp/time",n,"-",run_id,".RDS"))
    saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
    if(debug) { save(list = ls(), file = paste0("temp/all",n,"-",run_id,".rdata")); print(formulas) }
  }


  # End
  OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("Ending task. Selected: \n", formulas))
  stopCluster(cl)
  saveRDS(formulas, paste0("temp/featureselection_formulas",stamp,".RDS"))
  #saveRDS(formulas, paste0("featureselection_formulas.RDS"))
  #nazwa = paste0("temp/featureselection_data",stamp,".rda")
  #save(list = ls(), file = nazwa)
  #dev.off()
  print(formulas)
  setwd(oldwd)
  formulas}, timeout = timeout_sec, onTimeout = "silent")
  if (is.null(wynik_finalny)) {
    OmicSelector_log(logfile = "temp/featureselection.log",  message_to_log = paste0("STOPED BECAUSE OF ERROR OR TIMEOUT REACHED!! Debug result: ", print(wynik_finalny)))
    warnings()
    save(list = ls(), file = paste0("temp/timeoutdebug_",stamp,"-",m,".rdata"))
  }
  sink()
  #dev.off()
return(wynik_finalny)
}
