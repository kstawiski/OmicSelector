#' Mutual Information Feature Selection
#'
#' @description
#' Select features using mutual information to measure dependency between
#' features and class labels.
#'
#' @author OmicSelector Team

# =============================================================================
# METHOD IMPLEMENTATION
# =============================================================================

#' Mutual Information Feature Selection
#'
#' @description
#' Select features based on mutual information with the class variable.
#' Uses entropy-based measures to capture non-linear relationships.
#'
#' @param data Training data with class column as first column
#' @param config Configuration list (discretization method, etc.)
#' @param max_features Maximum number of features to select
#' @param use_smote Whether to apply SMOTE balancing
#' @param timeout_sec Timeout in seconds
#'
#' @return List with selected features and metadata
fs_method_mutual_information <- function(data, 
                                        config = list(), 
                                        max_features = 20, 
                                        use_smote = FALSE, 
                                        timeout_sec = 300) {
  
  # Check if required packages are available
  if (!requireNamespace("infotheo", quietly = TRUE) && 
      !requireNamespace("entropy", quietly = TRUE)) {
    
    # Fallback to basic implementation
    return(fs_method_mutual_information_basic(data, config, max_features, use_smote, timeout_sec))
  }
  
  # Apply SMOTE if requested
  if (use_smote) {
    data <- apply_smote_safely(data)
  }
  
  # Extract classes and features
  classes <- as.factor(data[, 1])
  features <- data[, -1, drop = FALSE]
  
  # Configuration parameters
  disc_method <- config$discretization %||% "equalwidth"
  n_bins <- config$n_bins %||% 5
  
  # Calculate mutual information
  mi_scores <- apply(features, 2, function(x) {
    tryCatch({
      
      if (requireNamespace("infotheo", quietly = TRUE)) {
        # Use infotheo package if available
        if (is.numeric(x)) {
          x_discrete <- infotheo::discretize(x, disc = disc_method, nbins = n_bins)
        } else {
          x_discrete <- as.factor(x)
        }
        
        mi_score <- infotheo::mutinformation(x_discrete, classes)
        return(mi_score)
        
      } else if (requireNamespace("entropy", quietly = TRUE)) {
        # Use entropy package as fallback
        if (is.numeric(x)) {
          x_discrete <- cut(x, breaks = n_bins, include.lowest = TRUE)
        } else {
          x_discrete <- as.factor(x)
        }
        
        # Calculate mutual information manually
        joint_entropy <- entropy::entropy(table(x_discrete, classes))
        x_entropy <- entropy::entropy(table(x_discrete))
        y_entropy <- entropy::entropy(table(classes))
        
        mi_score <- x_entropy + y_entropy - joint_entropy
        return(mi_score)
      }
      
      return(0)
      
    }, error = function(e) {
      return(0)
    })
  })
  
  # Handle edge cases
  mi_scores[is.na(mi_scores)] <- 0
  mi_scores[is.infinite(mi_scores)] <- 0
  mi_scores[mi_scores < 0] <- 0  # MI should be non-negative
  
  # Select top features
  if (length(mi_scores) == 0) {
    warning("No features could be evaluated for mutual information")
    return(NULL)
  }
  
  n_select <- min(max_features, length(mi_scores))
  selected_indices <- order(mi_scores, decreasing = TRUE)[seq_len(n_select)]
  selected_features <- names(features)[selected_indices]
  
  return(list(
    features = selected_features,
    scores = mi_scores[selected_indices],
    metadata = list(
      method = "Mutual Information",
      n_features_input = ncol(features),
      n_features_selected = length(selected_features),
      parameters_used = list(
        discretization = disc_method,
        n_bins = n_bins,
        max_features = max_features,
        use_smote = use_smote
      ),
      mi_scores = mi_scores[selected_indices]
    )
  ))
}

#' Basic Mutual Information Implementation
#'
#' @description
#' Fallback implementation that doesn't require external packages.
#'
#' @param data Training data
#' @param config Configuration
#' @param max_features Max features
#' @param use_smote SMOTE flag
#' @param timeout_sec Timeout
#'
#' @return Feature selection result
fs_method_mutual_information_basic <- function(data, config, max_features, use_smote, timeout_sec) {
  
  # Apply SMOTE if requested
  if (use_smote) {
    data <- apply_smote_safely(data)
  }
  
  # Extract classes and features
  classes <- as.factor(data[, 1])
  features <- data[, -1, drop = FALSE]
  
  # Configuration
  n_bins <- config$n_bins %||% 5
  
  # Basic entropy function
  entropy <- function(x) {
    probs <- table(x) / length(x)
    probs <- probs[probs > 0]  # Remove zero probabilities
    -sum(probs * log2(probs))
  }
  
  # Calculate mutual information scores
  mi_scores <- apply(features, 2, function(x) {
    tryCatch({
      # Discretize if numeric
      if (is.numeric(x)) {
        x_discrete <- cut(x, breaks = n_bins, include.lowest = TRUE, labels = FALSE)
      } else {
        x_discrete <- as.numeric(as.factor(x))
      }
      
      # Calculate entropies
      h_x <- entropy(x_discrete)
      h_y <- entropy(classes)
      h_xy <- entropy(paste(x_discrete, classes, sep = "_"))
      
      # Mutual information: I(X;Y) = H(X) + H(Y) - H(X,Y)
      mi <- h_x + h_y - h_xy
      
      return(max(0, mi))  # Ensure non-negative
      
    }, error = function(e) {
      return(0)
    })
  })
  
  # Handle edge cases
  mi_scores[is.na(mi_scores)] <- 0
  mi_scores[is.infinite(mi_scores)] <- 0
  
  # Select top features
  n_select <- min(max_features, length(mi_scores))
  selected_indices <- order(mi_scores, decreasing = TRUE)[seq_len(n_select)]
  selected_features <- names(features)[selected_indices]
  
  return(list(
    features = selected_features,
    scores = mi_scores[selected_indices],
    metadata = list(
      method = "Mutual Information (Basic)",
      n_features_input = ncol(features),
      n_features_selected = length(selected_features),
      parameters_used = list(
        n_bins = n_bins,
        max_features = max_features,
        use_smote = use_smote
      ),
      implementation = "basic_fallback"
    )
  ))
}

# =============================================================================
# METHOD REGISTRATION
# =============================================================================

# Register the method with the framework
tryCatch({
  register_fs_method(
    method_id = 102,
    method_name = "Mutual Information",
    method_function = fs_method_mutual_information,
    category = "filter",
    dependencies = c("infotheo", "entropy"),  # Optional packages
    description = "Select features using mutual information with class variable",
    parameters = list(
      discretization = "equalwidth",
      n_bins = 5
    ),
    complexity = "medium",
    supports_smote = TRUE,
    timeout_default = 180
  )
  
  cat("✓ Registered Mutual Information feature selection method (ID: 102)\n")
  
}, error = function(e) {
  warning("Failed to register Mutual Information method: ", e$message)
})
