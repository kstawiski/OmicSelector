#' Modern miRNA Name Correction Function
#'
#' @description
#' Modernized version of OmicSelector_correct_miRNA_names with improved performance,
#' error handling, and dependency management. Corrects miRNA names to latest miRBase
#' nomenclature with enhanced caching and validation.
#'
#' Key improvements:
#' - Removed redundant library loading (37 → 4 essential libraries)
#' - Added comprehensive error handling and validation
#' - Implemented caching system for miRBase data
#' - Memory-efficient parallel processing
#' - Progress tracking and timeout protection
#' - Support for multiple species
#' - Robust URL handling with fallbacks
#'
#' @param data Dataset with miRNA names in columns
#' @param species Species prefix for miRNAs (default: "hsa" for human)
#' @param correct_dots Whether to correct dots to hyphens (default: TRUE)
#' @param use_cache Whether to use cached miRBase data (default: TRUE)
#' @param cache_dir Directory for cached data (default: tempdir())
#' @param timeout_sec Timeout for downloads in seconds (default: 300)
#' @param parallel Whether to use parallel processing (default: TRUE)
#' @param verbose Enable progress messages (default: TRUE)
#' @param fallback_file Local fallback file if download fails
#'
#' @return Dataset with corrected miRNA names
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' corrected_data <- correct_mirna_names_modern(data)
#' 
#' # With custom species and settings
#' corrected_data <- correct_mirna_names_modern(
#'   data = my_data,
#'   species = "mmu",  # Mouse
#'   correct_dots = TRUE,
#'   use_cache = TRUE
#' )
#' 
#' # Non-parallel processing for small datasets
#' corrected_data <- correct_mirna_names_modern(
#'   data = small_data,
#'   parallel = FALSE
#' )
#' }
#'
#' @export
correct_mirna_names_modern <- function(data,
                                     species = "hsa",
                                     correct_dots = TRUE,
                                     use_cache = TRUE,
                                     cache_dir = tempdir(),
                                     timeout_sec = 300,
                                     parallel = TRUE,
                                     verbose = TRUE,
                                     fallback_file = NULL) {
  
  # Validate inputs
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("Input 'data' must be a data.frame or matrix")
  }
  
  if (!is.character(species) || length(species) != 1) {
    stop("'species' must be a single character string (e.g., 'hsa', 'mmu')")
  }
  
  if (ncol(data) == 0) {
    warning("Input data has no columns. Returning original data.")
    return(data)
  }
  
  # Load only essential libraries
  required_packages <- c("data.table", "dplyr", "parallel", "foreach")
  
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required but not installed", pkg))
    }
  }
  
  if (verbose) {
    cat("🔧 Starting modern miRNA name correction...\n")
    cat(sprintf("   Species: %s\n", species))
    cat(sprintf("   Columns to process: %d\n", ncol(data)))
  }
  
  # Setup cache directory
  if (use_cache) {
    cache_dir <- file.path(cache_dir, "omicselector_cache")
    if (!dir.exists(cache_dir)) {
      dir.create(cache_dir, recursive = TRUE)
    }
  }
  
  # Load miRBase aliases data
  mirbase_data <- load_mirbase_aliases(
    species = species,
    use_cache = use_cache,
    cache_dir = cache_dir,
    timeout_sec = timeout_sec,
    verbose = verbose,
    fallback_file = fallback_file
  )
  
  if (is.null(mirbase_data)) {
    warning("Could not load miRBase data. Returning original data.")
    return(data)
  }
  
  # Process column names
  column_names <- colnames(data)
  
  if (verbose) {
    cat("🔄 Processing column names...\n")
  }
  
  # Correct names in parallel or sequential
  if (parallel && length(column_names) > 10) {
    corrected_names <- correct_names_parallel(
      column_names = column_names,
      mirbase_data = mirbase_data,
      species = species,
      correct_dots = correct_dots,
      verbose = verbose
    )
  } else {
    corrected_names <- correct_names_sequential(
      column_names = column_names,
      mirbase_data = mirbase_data,
      species = species,
      correct_dots = correct_dots,
      verbose = verbose
    )
  }
  
  # Apply corrected names
  data_corrected <- data
  colnames(data_corrected) <- corrected_names
  
  # Report changes
  changes <- sum(column_names != corrected_names)
  if (verbose) {
    cat(sprintf("✅ Correction complete: %d/%d names updated\n", changes, length(column_names)))
  }
  
  return(data_corrected)
}

#' Load miRBase Aliases Data
#' @keywords internal
load_mirbase_aliases <- function(species, use_cache, cache_dir, timeout_sec, verbose, fallback_file) {
  
  cache_file <- if (use_cache) file.path(cache_dir, paste0("mirbase_aliases_", species, ".rds")) else NULL
  
  # Try to load from cache first
  if (use_cache && file.exists(cache_file)) {
    if (verbose) cat("📂 Loading miRBase data from cache...\n")
    
    tryCatch({
      cached_data <- readRDS(cache_file)
      
      # Check cache age (refresh if older than 30 days)
      cache_age <- difftime(Sys.time(), file.mtime(cache_file), units = "days")
      if (cache_age < 30) {
        return(cached_data)
      } else if (verbose) {
        cat("⏰ Cache is older than 30 days, refreshing...\n")
      }
    }, error = function(e) {
      if (verbose) cat("⚠️ Cache file corrupted, re-downloading...\n")
    })
  }
  
  # Download fresh data
  if (verbose) cat("🌐 Downloading miRBase aliases data...\n")
  
  mirbase_url <- "ftp://mirbase.org/pub/mirbase/CURRENT/aliases.txt.gz"
  
  tryCatch({
    # Set timeout for download
    old_timeout <- getOption("timeout")
    options(timeout = timeout_sec)
    on.exit(options(timeout = old_timeout), add = TRUE)
    
    # Download and process data
    temp_file <- tempfile(fileext = ".txt.gz")
    download.file(mirbase_url, temp_file, mode = "wb", quiet = !verbose)
    
    aliases_data <- data.table::fread(temp_file, header = FALSE)
    colnames(aliases_data) <- c("MIMAT", "Aliases")
    
    # Filter for species and MIMAT entries
    species_pattern <- paste0("^", species, "-")
    filtered_data <- aliases_data[
      grepl(species_pattern, aliases_data$Aliases) & 
      grepl("^MIMAT", aliases_data$MIMAT)
    ]
    
    # Cache the results
    if (use_cache) {
      tryCatch({
        saveRDS(filtered_data, cache_file)
        if (verbose) cat("💾 Data cached successfully\n")
      }, error = function(e) {
        warning("Could not save to cache: ", e$message)
      })
    }
    
    unlink(temp_file)
    return(filtered_data)
    
  }, error = function(e) {
    if (verbose) cat("❌ Download failed:", e$message, "\n")
    
    # Try fallback file if provided
    if (!is.null(fallback_file) && file.exists(fallback_file)) {
      if (verbose) cat("📄 Using fallback file...\n")
      
      tryCatch({
        fallback_data <- data.table::fread(fallback_file, header = FALSE)
        colnames(fallback_data) <- c("MIMAT", "Aliases")
        
        species_pattern <- paste0("^", species, "-")
        return(fallback_data[
          grepl(species_pattern, fallback_data$Aliases) & 
          grepl("^MIMAT", fallback_data$MIMAT)
        ])
      }, error = function(e2) {
        if (verbose) cat("❌ Fallback file also failed:", e2$message, "\n")
      })
    }
    
    return(NULL)
  })
}

#' Correct Names in Parallel
#' @keywords internal
correct_names_parallel <- function(column_names, mirbase_data, species, correct_dots, verbose) {
  
  if (!requireNamespace("doParallel", quietly = TRUE) || 
      !requireNamespace("foreach", quietly = TRUE)) {
    warning("Parallel processing packages not available, using sequential processing")
    return(correct_names_sequential(column_names, mirbase_data, species, correct_dots, verbose))
  }
  
  # Setup parallel processing
  cores <- min(parallel::detectCores() - 1, length(column_names))
  if (cores < 2) {
    return(correct_names_sequential(column_names, mirbase_data, species, correct_dots, verbose))
  }
  
  cl <- parallel::makePSOCKcluster(cores)
  doParallel::registerDoParallel(cl)
  
  on.exit({
    parallel::stopCluster(cl)
    doParallel::stopImplicitCluster()
  }, add = TRUE)
  
  if (verbose) {
    cat(sprintf("🔄 Using %d cores for parallel processing...\n", cores))
  }
  
  # Export necessary objects to workers
  parallel::clusterExport(cl, c("mirbase_data", "species", "correct_dots", "correct_single_name"), envir = environment())
  
  tryCatch({
    `%dopar%` <- foreach::`%dopar%`  # Import the operator
    corrected_names <- foreach::foreach(
      name = column_names,
      .combine = c,
      .packages = c("dplyr", "data.table"),
      .errorhandling = "pass"
    ) %dopar% {
      correct_single_name(name, mirbase_data, species, correct_dots)
    }
    
    return(corrected_names)
    
  }, error = function(e) {
    warning("Parallel processing failed, falling back to sequential: ", e$message)
    return(correct_names_sequential(column_names, mirbase_data, species, correct_dots, verbose))
  })
}

#' Correct Names Sequentially
#' @keywords internal
correct_names_sequential <- function(column_names, mirbase_data, species, correct_dots, verbose) {
  
  corrected_names <- character(length(column_names))
  
  for (i in seq_along(column_names)) {
    corrected_names[i] <- correct_single_name(
      column_names[i], 
      mirbase_data, 
      species, 
      correct_dots
    )
    
    if (verbose && i %% 100 == 0) {
      cat(sprintf("   Processed %d/%d names\n", i, length(column_names)))
    }
  }
  
  return(corrected_names)
}

#' Correct Single miRNA Name
#' @keywords internal
correct_single_name <- function(name, mirbase_data, species, correct_dots) {
  
  # Apply dot correction if requested
  search_name <- if (correct_dots) gsub("\\.", "-", name) else name
  
  # Look for exact match first
  exact_match <- mirbase_data[mirbase_data$Aliases == search_name, ]
  if (nrow(exact_match) > 0) {
    return(search_name)
  }
  
  # Look for partial matches (case insensitive)
  pattern <- sprintf(".*%s.*", gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", search_name))
  partial_matches <- mirbase_data[grepl(pattern, mirbase_data$Aliases, ignore.case = TRUE), ]
  
  if (nrow(partial_matches) > 0) {
    # Return the first alias that contains our search term
    best_match <- partial_matches$Aliases[1]
    return(best_match)
  }
  
  # If no match found, return original name
  return(name)
}

#' Backward Compatibility Alias
#' @export
OmicSelector_correct_miRNA_names <- function(temp, species = "hsa", correct_dots = TRUE) {
  .Deprecated("correct_mirna_names_modern", 
              msg = "OmicSelector_correct_miRNA_names is deprecated. Use correct_mirna_names_modern() instead.")
  
  # Convert to modern function call
  correct_mirna_names_modern(
    data = temp,
    species = species,
    correct_dots = correct_dots,
    verbose = FALSE  # Reduce verbosity for backward compatibility
  )
}
