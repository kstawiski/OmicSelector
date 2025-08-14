#' Modern Count Normalization for miRNA Expression Data
#'
#' @description
#' Converts raw miRNA counts to log10-transformed TPM (Transcripts Per Million) 
#' values with optional expression filtering. Designed specifically for miRNA 
#' data with improved error handling and modern R practices.
#'
#' @param count_data Matrix or data frame with miRNA counts (features in columns, samples in rows)
#' @param sample_metadata Data frame with sample metadata including Class column
#' @param sample_ids Character vector of sample identifiers (default: rownames of count_data)
#' @param apply_filter Logical indicating whether to apply expression filtering (default: TRUE)
#' @param min_counts Numeric minimum count threshold for filtering (default: 10)
#' @param min_samples_prop Numeric proportion of samples that must have min_counts (default: 0.5)
#' @param log_increment Numeric value added before log transformation (default: 0.001)
#' @param feature_prefix Character string prefix for valid features (default: "hsa")
#' @param validate_features Logical indicating whether to validate feature names (default: TRUE)
#' @param save_intermediate Logical indicating whether to save intermediate results (default: TRUE)
#' @param output_dir Character string for output directory (default: current directory)
#' @param verbose Logical indicating whether to print progress messages (default: TRUE)
#'
#' @return List containing:
#'   \item{normalized_data}{Matrix of log10-transformed TPM values}
#'   \item{tpm_data}{Matrix of TPM values before log transformation}
#'   \item{filtering_stats}{List with filtering statistics}
#'   \item{dgelist}{DGEList object for use with edgeR (if edgeR available)}
#'
#' @examples
#' \dontrun{
#' # Create sample data
#' counts <- matrix(rpois(1000, lambda = 50), nrow = 20, ncol = 50)
#' colnames(counts) <- paste0("hsa_miR_", 1:50)
#' rownames(counts) <- paste0("Sample_", 1:20)
#' 
#' metadata <- data.frame(
#'   Class = factor(rep(c("Case", "Control"), each = 10)),
#'   row.names = rownames(counts)
#' )
#' 
#' # Normalize counts
#' result <- normalize_mirna_counts(counts, metadata)
#' head(result$normalized_data[, 1:5])
#' }
#'
#' @export
normalize_mirna_counts <- function(count_data,
                                 sample_metadata = NULL,
                                 sample_ids = NULL,
                                 apply_filter = TRUE,
                                 min_counts = 10,
                                 min_samples_prop = 0.5,
                                 log_increment = 0.001,
                                 feature_prefix = "hsa",
                                 validate_features = TRUE,
                                 save_intermediate = TRUE,
                                 output_dir = ".",
                                 verbose = TRUE) {
  
  # Input validation
  if (!is.matrix(count_data) && !is.data.frame(count_data)) {
    stop("'count_data' must be a matrix or data frame", call. = FALSE)
  }
  
  if (is.data.frame(count_data)) {
    # Check if all columns are numeric
    numeric_cols <- sapply(count_data, is.numeric)
    if (!all(numeric_cols)) {
      if (verbose) message("Converting non-numeric columns to numeric")
      count_data[, !numeric_cols] <- lapply(count_data[, !numeric_cols], as.numeric)
    }
    count_data <- as.matrix(count_data)
  }
  
  # Validate feature names if requested
  if (validate_features && !is.null(feature_prefix)) {
    valid_features <- grepl(paste0("^", feature_prefix), colnames(count_data))
    if (!all(valid_features)) {
      invalid_count <- sum(!valid_features)
      if (verbose) {
        message("Removing ", invalid_count, " features that don't start with '", feature_prefix, "'")
      }
      count_data <- count_data[, valid_features, drop = FALSE]
    }
  }
  
  # Set sample IDs
  if (is.null(sample_ids)) {
    sample_ids <- rownames(count_data)
    if (is.null(sample_ids)) {
      sample_ids <- paste0("Sample_", seq_len(nrow(count_data)))
      rownames(count_data) <- sample_ids
    }
  }
  
  # Validate dimensions
  if (nrow(count_data) != length(sample_ids)) {
    stop("Number of rows in count_data must match length of sample_ids", call. = FALSE)
  }
  
  # Initialize filtering stats
  filtering_stats <- list(
    initial_features = ncol(count_data),
    initial_samples = nrow(count_data),
    features_after_filter = ncol(count_data),
    samples_after_filter = nrow(count_data),
    filter_applied = apply_filter
  )
  
  if (verbose) {
    message("Starting normalization of ", nrow(count_data), " samples and ", 
            ncol(count_data), " features")
  }
  
  # Apply expression filtering
  if (apply_filter) {
    min_samples <- ceiling(nrow(count_data) * min_samples_prop)
    
    if (verbose) {
      message("Applying filter: minimum ", min_counts, " counts in at least ", 
              min_samples, " samples (", round(min_samples_prop * 100, 1), "%)")
    }
    
    # Count samples meeting threshold for each feature
    feature_counts <- colSums(count_data >= min_counts)
    keep_features <- feature_counts >= min_samples
    
    if (sum(keep_features) == 0) {
      stop("No features pass the filtering criteria", call. = FALSE)
    }
    
    count_data <- count_data[, keep_features, drop = FALSE]
    filtering_stats$features_after_filter <- ncol(count_data)
    
    if (verbose) {
      message("Kept ", sum(keep_features), " features after filtering (",
              sum(!keep_features), " removed)")
    }
  }
  
  # Calculate library sizes (total counts per sample)
  lib_sizes <- rowSums(count_data)
  
  if (any(lib_sizes == 0)) {
    zero_lib_samples <- sum(lib_sizes == 0)
    warning("Found ", zero_lib_samples, " samples with zero library size", call. = FALSE)
    
    # Remove samples with zero library size
    keep_samples <- lib_sizes > 0
    count_data <- count_data[keep_samples, , drop = FALSE]
    lib_sizes <- lib_sizes[keep_samples]
    sample_ids <- sample_ids[keep_samples]
    
    filtering_stats$samples_after_filter <- nrow(count_data)
    
    if (verbose) {
      message("Removed ", zero_lib_samples, " samples with zero library size")
    }
  }
  
  # Calculate TPM values
  # TPM = (counts / feature_length) * 1e6 / (sum(counts / feature_length))
  # For miRNAs, assume uniform length (set to 1), so TPM = counts * 1e6 / lib_size
  
  if (verbose) message("Calculating TPM values...")
  
  tpm_data <- sweep(count_data, 1, lib_sizes, "/") * 1e6
  
  # Add increment and apply log10 transformation
  if (verbose) message("Applying log10 transformation...")
  
  log_tpm_data <- log10(tpm_data + log_increment)
  
  # Create DGEList object if edgeR is available
  dgelist <- NULL
  if (requireNamespace("edgeR", quietly = TRUE)) {
    if (verbose) message("Creating DGEList object...")
    
    # Prepare group information
    if (!is.null(sample_metadata) && "Class" %in% colnames(sample_metadata)) {
      groups <- sample_metadata$Class[match(sample_ids, rownames(sample_metadata))]
      if (any(is.na(groups))) {
        warning("Some samples not found in metadata; using default group", call. = FALSE)
        groups[is.na(groups)] <- "Unknown"
      }
    } else {
      groups <- rep("Group1", length(sample_ids))
    }
    
    dgelist <- edgeR::DGEList(
      counts = t(count_data),  # DGEList expects features in rows
      group = groups,
      samples = data.frame(
        sample_id = sample_ids,
        lib.size = lib_sizes,
        row.names = sample_ids
      )
    )
    
    # Calculate normalization factors
    dgelist <- edgeR::calcNormFactors(dgelist)
  }
  
  # Save intermediate results if requested
  if (save_intermediate) {
    if (verbose) message("Saving intermediate results to ", output_dir)
    
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Save TPM data
    saveRDS(tpm_data, file.path(output_dir, "TPM_data.rds"))
    saveRDS(log_tpm_data, file.path(output_dir, "log10_TPM_data.rds"))
    
    # Save DGEList if available
    if (!is.null(dgelist)) {
      saveRDS(dgelist, file.path(output_dir, "TPM_DGEList_filtered.rds"))
    }
    
    # Save filtering stats
    saveRDS(filtering_stats, file.path(output_dir, "normalization_stats.rds"))
  }
  
  if (verbose) {
    message("Normalization completed successfully!")
    message("Final dimensions: ", nrow(log_tpm_data), " samples x ", ncol(log_tpm_data), " features")
  }
  
  # Return results
  structure(list(
    normalized_data = log_tpm_data,
    tpm_data = tpm_data,
    filtering_stats = filtering_stats,
    dgelist = dgelist,
    sample_ids = sample_ids
  ), class = c("mirna_normalized", "list"))
}

#' @rdname normalize_mirna_counts  
#' @export
OmicSelector_counts_to_log10tpm <- function(danex, 
                                           metadane = NULL, 
                                           ids = NULL, 
                                           filtr = TRUE,
                                           filtr_minimalcounts = 10, 
                                           filtr_howmany = 1/2, 
                                           increment = 0.001) {
  
  # Backward compatibility wrapper
  .Deprecated("normalize_mirna_counts", 
              msg = "OmicSelector_counts_to_log10tpm is deprecated. Use normalize_mirna_counts instead.")
  
  # Map legacy parameters to new function
  if (is.null(ids) && !is.null(metadane)) {
    if ("ID" %in% colnames(metadane)) {
      ids <- metadane$ID
    }
  }
  
  result <- normalize_mirna_counts(
    count_data = danex,
    sample_metadata = metadane,
    sample_ids = ids,
    apply_filter = filtr,
    min_counts = filtr_minimalcounts,
    min_samples_prop = filtr_howmany,
    log_increment = increment,
    verbose = TRUE
  )
  
  # Return in legacy format (just the normalized data)
  return(result$normalized_data)
}

#' Print method for mirna_normalized objects
#' @export
print.mirna_normalized <- function(x, ...) {
  cat("miRNA Normalization Results\n")
  cat("===========================\n\n")
  
  cat("Dimensions:\n")
  cat("  Samples:", nrow(x$normalized_data), "\n")
  cat("  Features:", ncol(x$normalized_data), "\n\n")
  
  cat("Filtering Statistics:\n")
  cat("  Initial features:", x$filtering_stats$initial_features, "\n")
  cat("  Features after filter:", x$filtering_stats$features_after_filter, "\n")
  cat("  Filter applied:", x$filtering_stats$filter_applied, "\n\n")
  
  cat("Data Summary (log10 TPM):\n")
  print(summary(as.vector(x$normalized_data)))
  
  invisible(x)
}
