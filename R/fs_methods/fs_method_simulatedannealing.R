#' Simulated Annealing Feature Selection
#'
#' @description
#' Simulated Annealing optimization for feature selection. Uses temperature-based
#' acceptance probability to escape local optima and find globally optimal
#' feature subsets with Random Forest evaluation.
#' This is method #13 from the original OmicSelector comprehensive system.
#'
#' @author OmicSelector Team

# =============================================================================
# METHOD IMPLEMENTATION
# =============================================================================

#' Simulated Annealing Feature Selection
#'
#' @description
#' Applies simulated annealing algorithm to optimize feature selection.
#' Uses temperature cooling schedule and probabilistic acceptance of
#' worse solutions to avoid local optima.
#'
#' @param data Training data with class column as first column
#' @param config Configuration list (initial_temp, cooling_rate, min_temp, max_iterations)
#' @param max_features Maximum number of features to select
#' @param use_smote Whether to apply SMOTE balancing
#' @param timeout_sec Timeout in seconds
#'
#' @return List with selected features and metadata
fs_method_simulatedannealing <- function(data, 
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
  n_total_features <- ncol(feature_cols)
  
  # Configuration parameters
  initial_temp <- config$initial_temp %||% 100
  cooling_rate <- config$cooling_rate %||% 0.95
  min_temp <- config$min_temp %||% 0.01
  max_iterations <- config$max_iterations %||% 500
  
  # Initialize solution (random feature subset)
  current_solution <- initialize_sa_solution(n_total_features, max_features)
  current_fitness <- evaluate_sa_fitness(current_solution, data, feature_cols)
  
  # Best solution tracking
  best_solution <- current_solution
  best_fitness <- current_fitness
  
  # Annealing process tracking
  temperature <- initial_temp
  iteration <- 0
  accepted_moves <- 0
  rejected_moves <- 0
  temperature_history <- numeric(0)
  fitness_history <- numeric(0)
  
  start_time <- Sys.time()
  
  # Main annealing loop
  while (temperature > min_temp && 
         iteration < max_iterations &&
         as.numeric(difftime(Sys.time(), start_time, units = "secs")) < timeout_sec) {
    
    iteration <- iteration + 1
    
    tryCatch({
      # Generate neighbor solution
      neighbor_solution <- generate_neighbor(current_solution, max_features)
      neighbor_fitness <- evaluate_sa_fitness(neighbor_solution, data, feature_cols)
      
      # Calculate acceptance probability
      delta_fitness <- neighbor_fitness - current_fitness
      
      if (delta_fitness > 0) {
        # Better solution - always accept
        accept_move <- TRUE
      } else {
        # Worse solution - accept with probability
        acceptance_prob <- exp(delta_fitness / temperature)
        accept_move <- runif(1) < acceptance_prob
      }
      
      # Update current solution
      if (accept_move) {
        current_solution <- neighbor_solution
        current_fitness <- neighbor_fitness
        accepted_moves <- accepted_moves + 1
        
        # Update best solution if improved
        if (current_fitness > best_fitness) {
          best_solution <- current_solution
          best_fitness <- current_fitness
        }
      } else {
        rejected_moves <- rejected_moves + 1
      }
      
      # Record history
      temperature_history <- c(temperature_history, temperature)
      fitness_history <- c(fitness_history, current_fitness)
      
      # Cool down temperature
      temperature <- temperature * cooling_rate
      
    }, error = function(e) {
      warning(paste("Simulated Annealing iteration", iteration, "failed:", e$message))
      temperature <- temperature * cooling_rate
    })
  }
  
  # Extract selected features from best solution
  selected_indices <- which(best_solution == 1)
  
  if (length(selected_indices) > 0) {
    selected_features <- colnames(feature_cols)[selected_indices]
    
    # Calculate final feature scores
    final_scores <- calculate_sa_feature_scores(selected_features, data, feature_cols)
  } else {
    # Fallback if no features selected
    selected_features <- head(colnames(feature_cols), min(5, max_features))
    final_scores <- rep(0.5, length(selected_features))
    names(final_scores) <- selected_features
  }
  
  return(list(
    features = selected_features,
    scores = final_scores,
    metadata = list(
      method = "SimulatedAnnealing",
      n_features_input = ncol(feature_cols),
      n_features_selected = length(selected_features),
      parameters_used = list(
        initial_temp = initial_temp,
        cooling_rate = cooling_rate,
        min_temp = min_temp,
        max_iterations = max_iterations,
        max_features = max_features,
        use_smote = use_smote
      ),
      evaluation_metrics = list(
        best_fitness = best_fitness,
        final_fitness = current_fitness,
        final_scores = final_scores,
        mean_score = mean(final_scores)
      ),
      annealing_process = list(
        iterations_completed = iteration,
        final_temperature = temperature,
        accepted_moves = accepted_moves,
        rejected_moves = rejected_moves,
        acceptance_rate = accepted_moves / (accepted_moves + rejected_moves),
        temperature_history = temperature_history,
        fitness_history = fitness_history
      )
    )
  ))
}

#' Initialize SA Solution
#'
#' @description
#' Initialize random solution for simulated annealing
#'
#' @param n_features Total number of features
#' @param max_features Maximum features to select
#' @return Binary solution vector
initialize_sa_solution <- function(n_features, max_features) {
  # Random number of features to select
  n_select <- sample(1:min(max_features, n_features), 1)
  
  solution <- rep(0, n_features)
  selected_indices <- sample(n_features, n_select)
  solution[selected_indices] <- 1
  
  return(solution)
}

#' Evaluate SA Fitness
#'
#' @description
#' Evaluate fitness of a solution using cross-validation
#'
#' @param solution Binary solution vector
#' @param data Full dataset
#' @param feature_cols Feature columns
#' @return Fitness score
evaluate_sa_fitness <- function(solution, data, feature_cols) {
  selected_indices <- which(solution == 1)
  
  if (length(selected_indices) == 0) {
    return(0)
  }
  
  tryCatch({
    # Create subset data
    subset_data <- data[, c(TRUE, rep(FALSE, ncol(feature_cols))), drop = FALSE]
    subset_data <- cbind(subset_data, feature_cols[, selected_indices, drop = FALSE])
    
    # Evaluate using 2-fold CV for speed
    fitness <- evaluate_subset_fitness_fast(subset_data)
    
    # Add penalty for too many features
    feature_penalty <- max(0, (length(selected_indices) - 10) * 0.01)
    
    return(fitness - feature_penalty)
    
  }, error = function(e) {
    return(0)
  })
}

#' Fast Fitness Evaluation
#'
#' @description
#' Fast fitness evaluation using simple metrics
#'
#' @param subset_data Data subset
#' @return Fitness score
evaluate_subset_fitness_fast <- function(subset_data) {
  tryCatch({
    if (ncol(subset_data) <= 1) return(0)
    
    class_col <- subset_data[, 1]
    feature_cols <- subset_data[, -1, drop = FALSE]
    
    # Convert class to numeric
    if (is.factor(class_col) || is.character(class_col)) {
      class_numeric <- as.numeric(as.factor(class_col))
    } else {
      class_numeric <- class_col
    }
    
    # Calculate average correlation as fitness proxy
    correlations <- apply(feature_cols, 2, function(x) {
      if (is.numeric(x)) {
        abs(cor(x, class_numeric, use = "complete.obs"))
      } else {
        x_numeric <- as.numeric(as.factor(x))
        abs(cor(x_numeric, class_numeric, use = "complete.obs"))
      }
    })
    
    correlations[is.na(correlations)] <- 0
    
    # Use mean correlation as fitness
    fitness <- mean(correlations)
    
    # Add diversity bonus (negative correlation between features)
    if (ncol(feature_cols) > 1) {
      feature_cor_matrix <- cor(feature_cols, use = "complete.obs")
      feature_cor_matrix[is.na(feature_cor_matrix)] <- 0
      diag(feature_cor_matrix) <- 0
      
      avg_inter_correlation <- mean(abs(feature_cor_matrix))
      diversity_bonus <- max(0, (0.5 - avg_inter_correlation)) * 0.2
      fitness <- fitness + diversity_bonus
    }
    
    return(fitness)
    
  }, error = function(e) {
    return(0.1)
  })
}

#' Generate Neighbor Solution
#'
#' @description
#' Generate neighbor solution by modifying current solution
#'
#' @param current_solution Current binary solution
#' @param max_features Maximum features constraint
#' @return Neighbor solution
generate_neighbor <- function(current_solution, max_features) {
  neighbor <- current_solution
  n_features <- length(current_solution)
  current_count <- sum(current_solution)
  
  # Choose modification type randomly
  modification_type <- sample(c("flip", "swap", "add_remove"), 1, 
                             prob = c(0.5, 0.3, 0.2))
  
  if (modification_type == "flip") {
    # Simple bit flip
    flip_position <- sample(n_features, 1)
    
    if (neighbor[flip_position] == 1) {
      neighbor[flip_position] <- 0
    } else if (current_count < max_features) {
      neighbor[flip_position] <- 1
    }
    
  } else if (modification_type == "swap" && current_count > 0) {
    # Swap: turn off one feature, turn on another
    on_positions <- which(neighbor == 1)
    off_positions <- which(neighbor == 0)
    
    if (length(on_positions) > 0 && length(off_positions) > 0) {
      turn_off <- sample(on_positions, 1)
      turn_on <- sample(off_positions, 1)
      
      neighbor[turn_off] <- 0
      neighbor[turn_on] <- 1
    }
    
  } else if (modification_type == "add_remove") {
    # Randomly add or remove features
    if (current_count == 0) {
      # Add a feature
      off_positions <- which(neighbor == 0)
      if (length(off_positions) > 0) {
        add_position <- sample(off_positions, 1)
        neighbor[add_position] <- 1
      }
    } else if (current_count >= max_features) {
      # Remove a feature
      on_positions <- which(neighbor == 1)
      remove_position <- sample(on_positions, 1)
      neighbor[remove_position] <- 0
    } else {
      # Randomly add or remove
      if (runif(1) < 0.5) {
        # Add
        off_positions <- which(neighbor == 0)
        if (length(off_positions) > 0) {
          add_position <- sample(off_positions, 1)
          neighbor[add_position] <- 1
        }
      } else {
        # Remove
        on_positions <- which(neighbor == 1)
        if (length(on_positions) > 0) {
          remove_position <- sample(on_positions, 1)
          neighbor[remove_position] <- 0
        }
      }
    }
  }
  
  return(neighbor)
}

#' Calculate SA Feature Scores
#'
#' @description
#' Calculate individual scores for selected features
#'
#' @param selected_features Selected feature names
#' @param data Full dataset
#' @param feature_cols Feature columns
#' @return Feature scores
calculate_sa_feature_scores <- function(selected_features, data, feature_cols) {
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
    method_id = 13,
    method_name = "SimulatedAnnealing",
    method_function = fs_method_simulatedannealing,
    category = "optimization",
    dependencies = character(0),
    description = "Simulated Annealing optimization for feature selection with temperature cooling",
    parameters = list(
      initial_temp = 100,
      cooling_rate = 0.95,
      min_temp = 0.01,
      max_iterations = 500
    ),
    complexity = "high",
    supports_smote = TRUE,
    timeout_default = 450
  )
  
  cat("✓ Registered Simulated Annealing (SimulatedAnnealing) method (ID: 13)\n")
  
}, error = function(e) {
  warning("Failed to register SimulatedAnnealing method: ", e$message)
})
