#' Tabu Search Feature Selection
#'
#' @description
#' Tabu Search optimization for feature selection. Uses memory-based search
#' with tabu list to avoid cycling and explore feature space efficiently
#' while preventing revisiting recent solutions.
#' This is method #16 from the original OmicSelector comprehensive system.
#'
#' @author OmicSelector Team

# =============================================================================
# METHOD IMPLEMENTATION
# =============================================================================

#' Tabu Search Feature Selection
#'
#' @description
#' Implements Tabu Search algorithm for feature selection with adaptive
#' tabu list management and diversification strategies to find optimal
#' feature subsets while avoiding local optima.
#'
#' @param data Training data with class column as first column
#' @param config Configuration list (max_iterations, tabu_length, intensification_threshold)
#' @param max_features Maximum number of features to select
#' @param use_smote Whether to apply SMOTE balancing
#' @param timeout_sec Timeout in seconds
#'
#' @return List with selected features and metadata
fs_method_tabu_search <- function(data, 
                                 config = list(), 
                                 max_features = 20, 
                                 use_smote = FALSE, 
                                 timeout_sec = 300) {
  
  # Apply SMOTE if requested
  if (use_smote) {
    data <- apply_smote_safely(data)
  }
  
  # Extract features
  feature_cols <- data[, -1, drop = FALSE]
  n_features <- ncol(feature_cols)
  
  # Configuration parameters
  max_iterations <- config$max_iterations %||% 200
  tabu_length <- config$tabu_length %||% 20
  intensification_threshold <- config$intensification_threshold %||% 50
  diversification_threshold <- config$diversification_threshold %||% 100
  
  # Initialize solution
  current_solution <- initialize_tabu_solution(n_features, max_features)
  current_fitness <- evaluate_tabu_fitness(current_solution, data, feature_cols)
  
  # Best solution tracking
  best_solution <- current_solution
  best_fitness <- current_fitness
  
  # Tabu search components
  tabu_list <- list()
  tabu_times <- numeric(0)
  frequency_matrix <- rep(0, n_features)  # Feature usage frequency
  
  # Search history
  search_history <- list()
  stagnation_counter <- 0
  
  iteration <- 0
  start_time <- Sys.time()
  
  # Main tabu search loop
  while (iteration < max_iterations &&
         as.numeric(difftime(Sys.time(), start_time, units = "secs")) < timeout_sec) {
    
    iteration <- iteration + 1
    
    tryCatch({
      # Generate neighborhood solutions
      neighborhood <- generate_tabu_neighborhood(current_solution, max_features)
      
      # Evaluate neighborhood (excluding tabu moves)
      best_neighbor <- NULL
      best_neighbor_fitness <- -Inf
      best_neighbor_is_tabu <- FALSE
      
      for (neighbor in neighborhood) {
        # Check if move is tabu
        move_is_tabu <- is_tabu_move(current_solution, neighbor, tabu_list)
        neighbor_fitness <- evaluate_tabu_fitness(neighbor, data, feature_cols)
        
        # Apply aspiration criterion (accept tabu move if it's globally best)
        aspiration_criterion <- neighbor_fitness > best_fitness
        
        if (!move_is_tabu || aspiration_criterion) {
          if (neighbor_fitness > best_neighbor_fitness) {
            best_neighbor <- neighbor
            best_neighbor_fitness <- neighbor_fitness
            best_neighbor_is_tabu <- move_is_tabu
          }
        }
      }
      
      # Update solution
      if (!is.null(best_neighbor)) {
        # Add current move to tabu list
        move_description <- describe_move(current_solution, best_neighbor)
        tabu_list <- append(tabu_list, list(move_description), 0)
        tabu_times <- c(iteration, tabu_times)
        
        # Maintain tabu list size
        if (length(tabu_list) > tabu_length) {
          tabu_list <- tabu_list[1:tabu_length]
          tabu_times <- tabu_times[1:tabu_length]
        }
        
        # Update current solution
        current_solution <- best_neighbor
        current_fitness <- best_neighbor_fitness
        
        # Update frequency matrix
        selected_features <- which(current_solution == 1)
        frequency_matrix[selected_features] <- frequency_matrix[selected_features] + 1
        
        # Update best solution
        if (current_fitness > best_fitness) {
          best_solution <- current_solution
          best_fitness <- current_fitness
          stagnation_counter <- 0
        } else {
          stagnation_counter <- stagnation_counter + 1
        }
        
        # Record search history
        search_history[[iteration]] <- list(
          iteration = iteration,
          fitness = current_fitness,
          n_features = sum(current_solution),
          is_tabu = best_neighbor_is_tabu,
          stagnation = stagnation_counter
        )
        
        # Diversification strategy
        if (stagnation_counter >= diversification_threshold) {
          current_solution <- diversify_solution(frequency_matrix, n_features, max_features)
          current_fitness <- evaluate_tabu_fitness(current_solution, data, feature_cols)
          stagnation_counter <- 0
        }
      } else {
        # No valid neighbor found - force diversification
        current_solution <- diversify_solution(frequency_matrix, n_features, max_features)
        current_fitness <- evaluate_tabu_fitness(current_solution, data, feature_cols)
      }
      
    }, error = function(e) {
      warning(paste("Tabu Search iteration", iteration, "failed:", e$message))
    })
  }
  
  # Extract selected features
  selected_indices <- which(best_solution == 1)
  
  if (length(selected_indices) > 0) {
    selected_features <- colnames(feature_cols)[selected_indices]
    final_scores <- calculate_tabu_feature_scores(selected_features, data, feature_cols)
  } else {
    # Fallback
    selected_features <- head(colnames(feature_cols), min(5, max_features))
    final_scores <- rep(0.5, length(selected_features))
    names(final_scores) <- selected_features
  }
  
  return(list(
    features = selected_features,
    scores = final_scores,
    metadata = list(
      method = "TabuSearch",
      n_features_input = ncol(feature_cols),
      n_features_selected = length(selected_features),
      parameters_used = list(
        max_iterations = max_iterations,
        tabu_length = tabu_length,
        intensification_threshold = intensification_threshold,
        diversification_threshold = diversification_threshold,
        max_features = max_features,
        use_smote = use_smote
      ),
      evaluation_metrics = list(
        best_fitness = best_fitness,
        final_scores = final_scores,
        mean_score = mean(final_scores)
      ),
      search_process = list(
        iterations_completed = iteration,
        search_history = search_history,
        frequency_matrix = frequency_matrix,
        final_tabu_list_size = length(tabu_list),
        diversifications_performed = sum(sapply(search_history, function(x) x$stagnation == 0 && x$iteration > 1))
      )
    )
  ))
}

#' Initialize Tabu Solution
#'
#' @description
#' Initialize random solution for tabu search
#'
#' @param n_features Total number of features
#' @param max_features Maximum features to select
#' @return Binary solution vector
initialize_tabu_solution <- function(n_features, max_features) {
  n_select <- sample(floor(max_features * 0.3):max_features, 1)
  solution <- rep(0, n_features)
  selected_indices <- sample(n_features, n_select)
  solution[selected_indices] <- 1
  return(solution)
}

#' Evaluate Tabu Fitness
#'
#' @description
#' Fast fitness evaluation for tabu search
#'
#' @param solution Binary solution vector
#' @param data Full dataset
#' @param feature_cols Feature columns
#' @return Fitness score
evaluate_tabu_fitness <- function(solution, data, feature_cols) {
  selected_indices <- which(solution == 1)
  
  if (length(selected_indices) == 0) {
    return(0)
  }
  
  tryCatch({
    selected_features <- feature_cols[, selected_indices, drop = FALSE]
    class_col <- data[, 1]
    
    if (is.factor(class_col) || is.character(class_col)) {
      class_numeric <- as.numeric(as.factor(class_col))
    } else {
      class_numeric <- class_col
    }
    
    # Multi-criteria fitness
    correlation_fitness <- calculate_correlation_fitness(selected_features, class_numeric)
    diversity_fitness <- calculate_diversity_fitness(selected_features)
    size_penalty <- calculate_size_penalty(length(selected_indices), 10)
    
    total_fitness <- 0.6 * correlation_fitness + 0.3 * diversity_fitness - 0.1 * size_penalty
    
    return(max(0, total_fitness))
    
  }, error = function(e) {
    return(0)
  })
}

#' Calculate Correlation Fitness
#'
#' @description
#' Calculate fitness based on correlation with class
#'
#' @param features Selected feature matrix
#' @param class_numeric Numeric class vector
#' @return Correlation fitness score
calculate_correlation_fitness <- function(features, class_numeric) {
  correlations <- apply(features, 2, function(x) {
    if (is.numeric(x)) {
      abs(cor(x, class_numeric, use = "complete.obs"))
    } else {
      x_numeric <- as.numeric(as.factor(x))
      abs(cor(x_numeric, class_numeric, use = "complete.obs"))
    }
  })
  
  correlations[is.na(correlations)] <- 0
  return(mean(correlations))
}

#' Calculate Diversity Fitness
#'
#' @description
#' Calculate fitness based on feature diversity (low inter-correlation)
#'
#' @param features Selected feature matrix
#' @return Diversity fitness score
calculate_diversity_fitness <- function(features) {
  if (ncol(features) <= 1) {
    return(1)  # Maximum diversity for single feature
  }
  
  tryCatch({
    # Calculate feature correlation matrix
    feature_cor <- cor(features, use = "complete.obs")
    feature_cor[is.na(feature_cor)] <- 0
    diag(feature_cor) <- 0  # Remove self-correlations
    
    # Diversity is inverse of average absolute correlation
    avg_correlation <- mean(abs(feature_cor))
    diversity_score <- max(0, 1 - avg_correlation)
    
    return(diversity_score)
    
  }, error = function(e) {
    return(0.5)
  })
}

#' Calculate Size Penalty
#'
#' @description
#' Calculate penalty for feature set size
#'
#' @param n_features Number of selected features
#' @param target_size Target feature set size
#' @return Size penalty score
calculate_size_penalty <- function(n_features, target_size) {
  if (n_features <= target_size) {
    return(0)
  } else {
    return((n_features - target_size) / target_size)
  }
}

#' Generate Tabu Neighborhood
#'
#' @description
#' Generate neighborhood solutions by feature additions/removals
#'
#' @param current_solution Current binary solution
#' @param max_features Maximum features constraint
#' @return List of neighbor solutions
generate_tabu_neighborhood <- function(current_solution, max_features) {
  neighborhood <- list()
  current_count <- sum(current_solution)
  
  # Add neighbors (turn on features)
  if (current_count < max_features) {
    off_indices <- which(current_solution == 0)
    for (i in sample(off_indices, min(5, length(off_indices)))) {
      neighbor <- current_solution
      neighbor[i] <- 1
      neighborhood <- append(neighborhood, list(neighbor))
    }
  }
  
  # Remove neighbors (turn off features)
  if (current_count > 1) {
    on_indices <- which(current_solution == 1)
    for (i in sample(on_indices, min(5, length(on_indices)))) {
      neighbor <- current_solution
      neighbor[i] <- 0
      neighborhood <- append(neighborhood, list(neighbor))
    }
  }
  
  # Swap neighbors (simultaneous add/remove)
  on_indices <- which(current_solution == 1)
  off_indices <- which(current_solution == 0)
  
  if (length(on_indices) > 0 && length(off_indices) > 0) {
    n_swaps <- min(3, length(on_indices), length(off_indices))
    for (i in seq_len(n_swaps)) {
      neighbor <- current_solution
      off_idx <- sample(on_indices, 1)
      on_idx <- sample(off_indices, 1)
      neighbor[off_idx] <- 0
      neighbor[on_idx] <- 1
      neighborhood <- append(neighborhood, list(neighbor))
    }
  }
  
  return(neighborhood)
}

#' Check if Move is Tabu
#'
#' @description
#' Check if a move is in the tabu list
#'
#' @param current_solution Current solution
#' @param neighbor_solution Neighbor solution
#' @param tabu_list Current tabu list
#' @return TRUE if move is tabu
is_tabu_move <- function(current_solution, neighbor_solution, tabu_list) {
  if (length(tabu_list) == 0) {
    return(FALSE)
  }
  
  move_desc <- describe_move(current_solution, neighbor_solution)
  
  for (tabu_move in tabu_list) {
    if (identical(move_desc, tabu_move)) {
      return(TRUE)
    }
  }
  
  return(FALSE)
}

#' Describe Move
#'
#' @description
#' Create description of move between solutions
#'
#' @param from_solution Source solution
#' @param to_solution Target solution
#' @return Move description
describe_move <- function(from_solution, to_solution) {
  diff_indices <- which(from_solution != to_solution)
  
  added <- diff_indices[to_solution[diff_indices] == 1]
  removed <- diff_indices[to_solution[diff_indices] == 0]
  
  return(list(added = added, removed = removed))
}

#' Diversify Solution
#'
#' @description
#' Create diversified solution based on frequency matrix
#'
#' @param frequency_matrix Feature usage frequency
#' @param n_features Total number of features
#' @param max_features Maximum features to select
#' @return Diversified solution
diversify_solution <- function(frequency_matrix, n_features, max_features) {
  # Select features with low usage frequency
  n_select <- sample(floor(max_features * 0.5):max_features, 1)
  
  # Weight selection by inverse frequency
  weights <- 1 / (frequency_matrix + 1)
  selected_indices <- sample(n_features, n_select, prob = weights)
  
  solution <- rep(0, n_features)
  solution[selected_indices] <- 1
  
  return(solution)
}

#' Calculate Tabu Feature Scores
#'
#' @description
#' Calculate individual scores for selected features
#'
#' @param selected_features Selected feature names
#' @param data Full dataset
#' @param feature_cols Feature columns
#' @return Feature scores
calculate_tabu_feature_scores <- function(selected_features, data, feature_cols) {
  scores <- numeric(length(selected_features))
  names(scores) <- selected_features
  
  class_col <- data[, 1]
  if (is.factor(class_col) || is.character(class_col)) {
    class_numeric <- as.numeric(as.factor(class_col))
  } else {
    class_numeric <- class_col
  }
  
  for (i in seq_along(selected_features)) {
    feature_name <- selected_features[i]
    feature_data <- feature_cols[, feature_name]
    
    if (is.numeric(feature_data)) {
      scores[i] <- abs(cor(feature_data, class_numeric, use = "complete.obs"))
    } else {
      feature_numeric <- as.numeric(as.factor(feature_data))
      scores[i] <- abs(cor(feature_numeric, class_numeric, use = "complete.obs"))
    }
  }
  
  scores[is.na(scores)] <- 0.1
  return(scores)
}

# =============================================================================
# METHOD REGISTRATION
# =============================================================================

# Register the method with the framework
tryCatch({
  register_fs_method(
    method_id = 16,
    method_name = "TabuSearch",
    method_function = fs_method_tabu_search,
    category = "optimization",
    dependencies = character(0),
    description = "Tabu Search optimization with memory-based exploration and diversification",
    parameters = list(
      max_iterations = 200,
      tabu_length = 20,
      intensification_threshold = 50,
      diversification_threshold = 100
    ),
    complexity = "high",
    supports_smote = TRUE,
    timeout_default = 450
  )
  
  cat("✓ Registered Tabu Search (TabuSearch) method (ID: 16)\n")
  
}, error = function(e) {
  warning("Failed to register TabuSearch method: ", e$message)
})
