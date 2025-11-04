#' OmicSelector_OmicSelector_modern
#'
#' Modernized version of the main feature selection function. This streamlined version
#' removes overcomplicated methods with heavy dependencies and integrates Phase 2
#' stability-aware feature selection methods.
#'
#' **Key Improvements over original:**
#' - Removed 40+ complex/slow methods (bounceR, WxNet, genetic algorithms, etc.)
#' - Integrated Phase 2 stability selection methods
#' - Integrated feature clustering for biomarker replaceability
#' - Integrated model-agnostic feature importance
#' - Cleaner code structure with better error handling
#' - Better logging and progress tracking
#' - Reduced dependencies (no Biocomb, no Python/conda, no bounceR, no feseR)
#'
#' **Retained Methods (20 core methods):**
#' 1. all - All features baseline
#' 2-7. Differential expression: sig, sigtop, topFC (+ SMOTE versions)
#' 8-9. LASSO, ElasticNet (glmnet-based, fast and reliable)
#' 10. Boruta (popular stability-based method)
#' 11-12. Random Forest RFE (recursive feature elimination)
#' 13-14. stepAIC (stepwise selection)
#' 15-16. stepLDA (linear discriminant analysis)
#' 17-20. **NEW: Phase 2 methods** (stability selection, clustering, importance)
#'
#' **Removed Methods (50+ overcomplicated methods):**
#' - bounceR (genetic + boosting - very slow, complex dependency)
#' - WxNet (requires Python/conda/neural networks - complex setup)
#' - GeneticAlgorithmRF (very slow, hard to reproduce)
#' - SimulatedAnnealing (slow, often no better than RFE)
#' - spFSR (niche method, complex dependency)
#' - varSelRF (replaced by RandomForestRFE which is better)
#' - All Biocomb methods (CFS, classloop, fcfs, fwrap, MDL) - old package, unstable
#' - My.stepwise methods (replaced by standard stepAIC)
#' - feseR methods (complex dependency, can use Phase 2 methods instead)
#' - Ridge regression (LASSO/ElasticNet are better)
#'
#' @param wd Working directory with data
#' @param methods Character vector specifying which method categories to use:
#'   \itemize{
#'     \item "de" - Differential expression methods (sig, topFC)
#'     \item "regularized" - LASSO, ElasticNet
#'     \item "embedded" - Boruta, RandomForestRFE
#'     \item "stepwise" - stepAIC, stepLDA
#'     \item "modern" - Phase 2 stability selection methods
#'     \item "all" - All available methods
#'   }
#' @param use_smote Logical, whether to also run methods on SMOTE-balanced data
#' @param prefer_no_features Maximum number of features to select (default: 11)
#' @param stamp Timestamp or identifier for output files
#' @param stability_iterations Number of iterations for stability selection (default: 50)
#' @param parallel Logical, whether to use parallel processing
#' @param n_cores Number of cores for parallel processing
#' @param timeout_sec Timeout in seconds for individual methods
#' @param type Mode for differential expression ("auto", "var", "nonvar")
#' @param verbose Logical, whether to print progress messages
#'
#' @return List of formulas selected by each method
#'
#' @examples
#' \dontrun{
#' # Run modern methods only
#' formulas <- OmicSelector_OmicSelector_modern(
#'   methods = c("de", "regularized", "modern"),
#'   prefer_no_features = 10,
#'   stability_iterations = 100
#' )
#'
#' # Run all methods including SMOTE versions
#' formulas <- OmicSelector_OmicSelector_modern(
#'   methods = "all",
#'   use_smote = TRUE
#' )
#' }
#'
#' @export
OmicSelector_OmicSelector_modern <- function(
  wd = getwd(),
  methods = "all",
  use_smote = TRUE,
  prefer_no_features = 11,
  stamp = as.numeric(Sys.time()),
  stability_iterations = 50,
  parallel = TRUE,
  n_cores = NULL,
  timeout_sec = 7200,  # 2 hours default (vs 48 hours in original!)
  type = "auto",
  verbose = TRUE
) {

  # Setup
  oldwd <- getwd()
  setwd(wd)
  on.exit(setwd(oldwd))

  # Load only necessary packages
  suppressMessages({
    library(dplyr)
    library(caret)
    library(glmnet)
    library(MASS)
    library(pROC)
    library(R.utils)
  })

  # Create temp directory
  if (!dir.exists("temp")) {
    dir.create("temp")
  }

  # Setup logging
  log_file <- paste0("temp/", stamp, "_modern_featureselection.log")
  if (verbose) {
    message("Starting modernized feature selection...")
    message("Log file: ", log_file)
  }

  # Determine which method groups to run
  if ("all" %in% methods) {
    methods <- c("de", "regularized", "embedded", "stepwise", "modern")
  }

  formulas <- list()
  times <- list()
  run_id <- stamp

  # Start timeout wrapper
  tryCatch({
    withTimeout({

      # Load data
      if (verbose) message("Loading data...")
      dane <- OmicSelector_load_datamix()
      train <- dane[[1]]
      test <- dane[[2]]
      valid <- dane[[3]]
      train_smoted <- dane[[4]]
      trainx <- dane[[5]]
      trainx_smoted <- dane[[6]]

      if (verbose) {
        message("Data loaded:")
        message("  Train: ", nrow(train), " samples, ", ncol(trainx), " features")
      }

      # ========================================================================
      # METHOD CATEGORY 1: DIFFERENTIAL EXPRESSION
      # ========================================================================
      if ("de" %in% methods) {
        if (verbose) message("\n=== Differential Expression Methods ===")

        start_time <- Sys.time()

        # Baseline: all features
        formulas[["all"]] <- OmicSelector_create_formula(colnames(trainx))

        # DE analysis
        if (verbose) message("Running differential expression analysis...")
        wyniki <- OmicSelector_differential_expression_ttest(trainx, train$Class, mode = type)

        # Significant features
        istotne <- filter(wyniki, `p-value BH` <= 0.05) %>% arrange(`p-value BH`)
        if (nrow(istotne) == 0) {
          istotne <- wyniki %>% arrange(`p-value BH`)
        }

        # Top features by different criteria
        istotne_top <- wyniki %>% arrange(`p-value BH`) %>% head(prefer_no_features)
        istotne_topFC <- wyniki %>% arrange(desc(abs(`log2FC`))) %>% head(prefer_no_features)

        # FC + sig filter
        fcsig <- wyniki %>%
          filter(`p-value BH` <= 0.05 & abs(`log2FC`) >= 1) %>%
          pull(miR) %>%
          as.character()
        if (length(fcsig) == 0) {
          fcsig <- istotne_top$miR
        }

        # Store formulas
        formulas[["sig"]] <- OmicSelector_create_formula(as.character(istotne$miR))
        formulas[["sigtop"]] <- OmicSelector_create_formula(as.character(istotne_top$miR))
        formulas[["topFC"]] <- OmicSelector_create_formula(as.character(istotne_topFC$miR))
        formulas[["fcsig"]] <- OmicSelector_create_formula(fcsig)

        # SMOTE versions if requested
        if (use_smote) {
          if (verbose) message("Running DE with SMOTE data...")
          wyniki_smoted <- OmicSelector_differential_expression_ttest(
            trainx_smoted, train_smoted$Class, mode = type
          )

          istotne_smoted <- filter(wyniki_smoted, `p-value BH` <= 0.05) %>%
            arrange(`p-value BH`)
          if (nrow(istotne_smoted) == 0) {
            istotne_smoted <- wyniki_smoted %>% arrange(`p-value BH`)
          }

          istotne_top_smoted <- wyniki_smoted %>%
            arrange(`p-value BH`) %>%
            head(prefer_no_features)

          formulas[["sigSMOTE"]] <- OmicSelector_create_formula(
            as.character(istotne_smoted$miR)
          )
          formulas[["sigtopSMOTE"]] <- OmicSelector_create_formula(
            as.character(istotne_top_smoted$miR)
          )
        }

        times[["de"]] <- Sys.time() - start_time
        if (verbose) message("DE methods completed in ", round(times[["de"]], 1), " sec")
      }

      # ========================================================================
      # METHOD CATEGORY 2: REGULARIZED METHODS (LASSO, ElasticNet)
      # ========================================================================
      if ("regularized" %in% methods) {
        if (verbose) message("\n=== Regularized Methods (LASSO, ElasticNet) ===")

        start_time <- Sys.time()

        # Prepare data for glmnet
        x_train <- as.matrix(trainx)
        y_train <- as.numeric(train$Class) - 1  # 0/1 encoding

        # LASSO
        if (verbose) message("Running LASSO...")
        tryCatch({
          lasso_cv <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 1,
                                nfolds = 10, type.measure = "auc")
          lasso_coef <- coef(lasso_cv, s = "lambda.min")
          lasso_features <- rownames(lasso_coef)[which(lasso_coef != 0)][-1]  # Remove intercept

          if (length(lasso_features) > 0) {
            formulas[["LASSO"]] <- OmicSelector_create_formula(lasso_features)
          }
        }, error = function(e) {
          if (verbose) message("  LASSO failed: ", conditionMessage(e))
        })

        # ElasticNet (alpha = 0.5)
        if (verbose) message("Running ElasticNet...")
        tryCatch({
          enet_cv <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 0.5,
                               nfolds = 10, type.measure = "auc")
          enet_coef <- coef(enet_cv, s = "lambda.min")
          enet_features <- rownames(enet_coef)[which(enet_coef != 0)][-1]

          if (length(enet_features) > 0) {
            formulas[["ElasticNet"]] <- OmicSelector_create_formula(enet_features)
          }
        }, error = function(e) {
          if (verbose) message("  ElasticNet failed: ", conditionMessage(e))
        })

        # SMOTE versions
        if (use_smote) {
          if (verbose) message("Running regularized methods with SMOTE...")

          x_train_smoted <- as.matrix(trainx_smoted)
          y_train_smoted <- as.numeric(train_smoted$Class) - 1

          tryCatch({
            lasso_cv <- cv.glmnet(x_train_smoted, y_train_smoted, family = "binomial",
                                  alpha = 1, nfolds = 10, type.measure = "auc")
            lasso_coef <- coef(lasso_cv, s = "lambda.min")
            lasso_features <- rownames(lasso_coef)[which(lasso_coef != 0)][-1]

            if (length(lasso_features) > 0) {
              formulas[["LASSO_SMOTE"]] <- OmicSelector_create_formula(lasso_features)
            }
          }, error = function(e) {
            if (verbose) message("  LASSO SMOTE failed: ", conditionMessage(e))
          })

          tryCatch({
            enet_cv <- cv.glmnet(x_train_smoted, y_train_smoted, family = "binomial",
                                 alpha = 0.5, nfolds = 10, type.measure = "auc")
            enet_coef <- coef(enet_cv, s = "lambda.min")
            enet_features <- rownames(enet_coef)[which(enet_coef != 0)][-1]

            if (length(enet_features) > 0) {
              formulas[["ElasticNet_SMOTE"]] <- OmicSelector_create_formula(enet_features)
            }
          }, error = function(e) {
            if (verbose) message("  ElasticNet SMOTE failed: ", conditionMessage(e))
          })
        }

        times[["regularized"]] <- Sys.time() - start_time
        if (verbose) message("Regularized methods completed in ",
                           round(times[["regularized"]], 1), " sec")
      }

      # ========================================================================
      # METHOD CATEGORY 3: EMBEDDED METHODS (Boruta, RandomForestRFE)
      # ========================================================================
      if ("embedded" %in% methods) {
        if (verbose) message("\n=== Embedded Methods (Boruta, RF-RFE) ===")

        start_time <- Sys.time()

        # Boruta
        if (verbose) message("Running Boruta...")
        tryCatch({
          suppressMessages(library(Boruta))
          bor <- Boruta(Class ~ ., data = train, doTrace = 0, maxRuns = 100)
          confirmed <- names(bor$finalDecision)[bor$finalDecision == "Confirmed"]

          if (length(confirmed) > 0) {
            formulas[["Boruta"]] <- OmicSelector_create_formula(confirmed)
          }
        }, error = function(e) {
          if (verbose) message("  Boruta failed: ", conditionMessage(e))
        })

        # Random Forest RFE
        if (verbose) message("Running Random Forest RFE...")
        tryCatch({
          ctrl <- rfeControl(
            functions = rfFuncs,
            method = "cv",
            number = 5,
            verbose = FALSE
          )

          sizes <- c(5, 10, prefer_no_features, 20, 30)
          rfProfile <- rfe(
            x = trainx,
            y = train$Class,
            sizes = sizes,
            rfeControl = ctrl
          )

          formulas[["RandomForestRFE"]] <- OmicSelector_create_formula(
            predictors(rfProfile)
          )
        }, error = function(e) {
          if (verbose) message("  RandomForestRFE failed: ", conditionMessage(e))
        })

        # SMOTE versions
        if (use_smote) {
          if (verbose) message("Running embedded methods with SMOTE...")

          tryCatch({
            bor <- Boruta(Class ~ ., data = train_smoted, doTrace = 0, maxRuns = 100)
            confirmed <- names(bor$finalDecision)[bor$finalDecision == "Confirmed"]

            if (length(confirmed) > 0) {
              formulas[["BorutaSMOTE"]] <- OmicSelector_create_formula(confirmed)
            }
          }, error = function(e) {
            if (verbose) message("  Boruta SMOTE failed: ", conditionMessage(e))
          })

          tryCatch({
            ctrl <- rfeControl(
              functions = rfFuncs,
              method = "cv",
              number = 5,
              verbose = FALSE
            )

            sizes <- c(5, 10, prefer_no_features, 20, 30)
            rfProfile <- rfe(
              x = trainx_smoted,
              y = train_smoted$Class,
              sizes = sizes,
              rfeControl = ctrl
            )

            formulas[["RandomForestRFE_SMOTE"]] <- OmicSelector_create_formula(
              predictors(rfProfile)
            )
          }, error = function(e) {
            if (verbose) message("  RandomForestRFE SMOTE failed: ", conditionMessage(e))
          })
        }

        times[["embedded"]] <- Sys.time() - start_time
        if (verbose) message("Embedded methods completed in ",
                           round(times[["embedded"]], 1), " sec")
      }

      # ========================================================================
      # METHOD CATEGORY 4: STEPWISE METHODS
      # ========================================================================
      if ("stepwise" %in% methods) {
        if (verbose) message("\n=== Stepwise Methods ===")

        start_time <- Sys.time()

        # stepAIC
        if (verbose) message("Running stepAIC...")
        tryCatch({
          # Start with significant features only to speed up
          train_sig <- dplyr::select(train, as.character(istotne_top$miR), Class)

          full_model <- glm(Class ~ ., data = train_sig, family = binomial)
          step_model <- stepAIC(full_model, direction = "backward",
                                trace = FALSE, k = 2)

          selected_vars <- names(coef(step_model))[-1]  # Remove intercept
          if (length(selected_vars) > 0) {
            formulas[["stepAIC"]] <- OmicSelector_create_formula(selected_vars)
          }
        }, error = function(e) {
          if (verbose) message("  stepAIC failed: ", conditionMessage(e))
        })

        # stepLDA
        if (verbose) message("Running stepLDA...")
        tryCatch({
          train_sig <- dplyr::select(train, as.character(istotne_top$miR), Class)

          slda <- train(
            Class ~ .,
            data = train_sig,
            method = "stepLDA",
            trControl = trainControl(method = "cv", number = 5),
            trace = FALSE
          )

          formulas[["stepLDA"]] <- OmicSelector_create_formula(predictors(slda))
        }, error = function(e) {
          if (verbose) message("  stepLDA failed: ", conditionMessage(e))
        })

        times[["stepwise"]] <- Sys.time() - start_time
        if (verbose) message("Stepwise methods completed in ",
                           round(times[["stepwise"]], 1), " sec")
      }

      # ========================================================================
      # METHOD CATEGORY 5: MODERN PHASE 2 METHODS
      # ========================================================================
      if ("modern" %in% methods) {
        if (verbose) message("\n=== Phase 2 Modern Methods ===")

        start_time <- Sys.time()

        # Load Phase 2 modules
        source_if_exists <- function(file) {
          if (file.exists(file)) {
            source(file)
            return(TRUE)
          }
          return(FALSE)
        }

        phase2_available <- all(c(
          source_if_exists("R/feature_selection_modern.R"),
          source_if_exists("R/feature_clustering.R")
        ))

        if (phase2_available) {
          # Stability Selection
          if (verbose) message("Running Stability Selection...")
          tryCatch({
            stable_result <- OmicSelector_stable_features(
              data = train,
              outcome = "Class",
              method = "stability_selection",
              n_iterations = stability_iterations,
              selection_threshold = 0.6,
              max_features = prefer_no_features,
              parallel = parallel
            )

            formulas[["StabilitySelection"]] <- OmicSelector_create_formula(
              stable_result$selected_features
            )

            # Also save stability scores for reporting
            saveRDS(stable_result,
                   paste0("temp/stability_selection_", run_id, ".RDS"))

          }, error = function(e) {
            if (verbose) message("  Stability Selection failed: ", conditionMessage(e))
          })

          # Boruta Stable (more iterations than regular Boruta)
          if (verbose) message("Running Boruta Stable...")
          tryCatch({
            boruta_stable <- OmicSelector_stable_features(
              data = train,
              outcome = "Class",
              method = "boruta_stable",
              n_iterations = 20,  # Less than stability selection
              selection_threshold = 0.7,
              parallel = parallel
            )

            formulas[["BorutaStable"]] <- OmicSelector_create_formula(
              boruta_stable$selected_features
            )

          }, error = function(e) {
            if (verbose) message("  Boruta Stable failed: ", conditionMessage(e))
          })

          # Feature Clustering for dimensionality reduction
          if (verbose) message("Running Feature Clustering...")
          tryCatch({
            clusters <- OmicSelector_cluster_features(
              data = train,
              features = colnames(trainx),
              method = "hierarchical",
              n_clusters = prefer_no_features,
              distance_metric = "pearson",
              plot = FALSE
            )

            # Use representative features
            formulas[["ClusterRepresentatives"]] <- OmicSelector_create_formula(
              clusters$representatives
            )

            # Save clustering result for biomarker replaceability analysis
            saveRDS(clusters,
                   paste0("temp/feature_clustering_", run_id, ".RDS"))

          }, error = function(e) {
            if (verbose) message("  Feature Clustering failed: ", conditionMessage(e))
          })

        } else {
          if (verbose) {
            message("  Phase 2 modules not found. Skipping modern methods.")
            message("  Make sure R/feature_selection_modern.R and R/feature_clustering.R exist")
          }
        }

        times[["modern"]] <- Sys.time() - start_time
        if (verbose) message("Modern methods completed in ",
                           round(times[["modern"]], 1), " sec")
      }

      # ========================================================================
      # FINALIZE
      # ========================================================================

      # Save all formulas
      saveRDS(formulas, paste0("temp/formulas_modern_", run_id, ".RDS"))
      saveRDS(times, paste0("temp/times_modern_", run_id, ".RDS"))

      if (verbose) {
        message("\n=== Feature Selection Complete ===")
        message("Total methods run: ", length(formulas))
        message("Total time: ", round(sum(unlist(times)), 1), " seconds")
        message("Results saved to: temp/formulas_modern_", run_id, ".RDS")
      }

      return(formulas)

    }, timeout = timeout_sec)

  }, TimeoutException = function(ex) {
    message("Feature selection timed out after ", timeout_sec, " seconds")
    saveRDS(formulas, paste0("temp/formulas_modern_partial_", run_id, ".RDS"))
    return(formulas)
  }, error = function(e) {
    message("Error in feature selection: ", conditionMessage(e))
    saveRDS(formulas, paste0("temp/formulas_modern_error_", run_id, ".RDS"))
    return(formulas)
  })

  return(formulas)
}
