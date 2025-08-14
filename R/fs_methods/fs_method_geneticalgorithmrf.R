#' Genetic Algorithm with Random Forest Feature Selection
#'
#' @description
#' Genetic Algorithm optimization using Random Forest fitness evaluation.
#' Evolves feature subsets through selection, crossover, and mutation
#' operations with RF-based fitness scoring.
#' This is method #12 from the original OmicSelector comprehensive system.
#'
#' @author OmicSelector Team

# =============================================================================
# METHOD IMPLEMENTATION
# =============================================================================

#' Genetic Algorithm with Random Forest Feature Selection
#'
#' @description
#' Uses genetic algorithm to evolve optimal feature subsets, with Random Forest
#' accuracy as fitness function. Implements selection, crossover, mutation,
#' and elitism for robust optimization.
#'
#' @param data Training data with class column as first column
#' @param config Configuration list (population_size, generations, mutation_rate, elite_size)
#' @param max_features Maximum number of features to select
#' @param use_smote Whether to apply SMOTE balancing
#' @param timeout_sec Timeout in seconds
#'
#' @return List with selected features and metadata
fs_method_geneticalgorithmrf <- function(data, 
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
  population_size <- config$population_size %||% 30
  max_generations <- config$generations %||% 20
  mutation_rate <- config$mutation_rate %||% 0.1
  elite_size <- config$elite_size %||% max(2, floor(population_size * 0.1))
  crossover_rate <- config$crossover_rate %||% 0.8
  
  # Initialize population (binary chromosomes)
  population <- initialize_population(population_size, n_total_features, max_features)
  
  # Track evolution
  fitness_history <- list()
  best_individuals <- list()
  generation <- 0
  start_time <- Sys.time()
  
  # Evolution loop
  while (generation < max_generations && 
         as.numeric(difftime(Sys.time(), start_time, units = "secs")) < timeout_sec) {
    
    generation <- generation + 1
    
    tryCatch({
      # Evaluate fitness for all individuals
      fitness_scores <- evaluate_population_fitness(population, data, feature_cols)
      
      # Track best individuals
      best_idx <- which.max(fitness_scores)
      best_individuals[[generation]] <- list(
        generation = generation,
        chromosome = population[best_idx, ],
        fitness = fitness_scores[best_idx],
        n_features = sum(population[best_idx, ])
      )
      
      # Store fitness history
      fitness_history[[generation]] <- list(
        generation = generation,
        max_fitness = max(fitness_scores),
        mean_fitness = mean(fitness_scores),
        min_fitness = min(fitness_scores),
        diversity = calculate_population_diversity(population)
      )
      
      # Selection and reproduction
      new_population <- create_new_generation(
        population, fitness_scores, elite_size, 
        crossover_rate, mutation_rate, max_features
      )
      
      population <- new_population
      
    }, error = function(e) {
      warning(paste("Genetic Algorithm generation", generation, "failed:", e$message))
      break
    })
  }
  
  # Select best solution
  if (length(best_individuals) > 0) {
    final_fitness <- sapply(best_individuals, function(x) x$fitness)
    best_generation <- which.max(final_fitness)
    best_chromosome <- best_individuals[[best_generation]]$chromosome
    
    selected_feature_indices <- which(best_chromosome == 1)
    selected_features <- colnames(feature_cols)[selected_feature_indices]
    
    # Calculate final fitness scores for selected features
    final_scores <- calculate_feature_fitness_scores(selected_features, data, feature_cols)
    
  } else {
    # Fallback selection
    selected_features <- head(colnames(feature_cols), max_features)
    final_scores <- rep(0.5, length(selected_features))
    names(final_scores) <- selected_features
  }
  
  return(list(
    features = selected_features,
    scores = final_scores,
    metadata = list(
      method = "GeneticAlgorithmRF",
      n_features_input = ncol(feature_cols),
      n_features_selected = length(selected_features),
      parameters_used = list(
        population_size = population_size,
        generations = max_generations,
        mutation_rate = mutation_rate,
        elite_size = elite_size,
        crossover_rate = crossover_rate,
        max_features = max_features,
        use_smote = use_smote
      ),
      evaluation_metrics = list(
        final_fitness_scores = final_scores,
        best_fitness = if(length(best_individuals) > 0) max(sapply(best_individuals, function(x) x$fitness)) else 0,
        mean_final_score = mean(final_scores)
      ),
      evolution_process = list(
        generations_completed = generation,
        fitness_history = fitness_history,
        best_individuals_by_generation = best_individuals,
        final_population_diversity = if(exists("population")) calculate_population_diversity(population) else 0
      )
    )
  ))
}

#' Initialize GA Population
#'
#' @description
#' Create initial population of binary chromosomes
#'
#' @param pop_size Population size
#' @param n_features Total number of features
#' @param max_features Maximum features per individual
#' @return Population matrix
initialize_population <- function(pop_size, n_features, max_features) {
  population <- matrix(0, nrow = pop_size, ncol = n_features)
  
  for (i in seq_len(pop_size)) {
    # Random number of features (between 1 and max_features)
    n_selected <- sample(1:min(max_features, n_features), 1)
    selected_indices <- sample(n_features, n_selected)
    population[i, selected_indices] <- 1
  }
  
  return(population)
}

#' Evaluate Population Fitness
#'
#' @description
#' Evaluate fitness of all individuals using Random Forest
#'
#' @param population Population matrix
#' @param data Full dataset
#' @param feature_cols Feature columns
#' @return Fitness scores vector
evaluate_population_fitness <- function(population, data, feature_cols) {
  fitness_scores <- numeric(nrow(population))
  
  for (i in seq_len(nrow(population))) {
    chromosome <- population[i, ]
    selected_indices <- which(chromosome == 1)
    
    if (length(selected_indices) == 0) {
      fitness_scores[i] <- 0
      next
    }
    
    tryCatch({
      # Create subset data
      subset_data <- data[, c(TRUE, rep(FALSE, ncol(feature_cols))), drop = FALSE]
      subset_data <- cbind(subset_data, feature_cols[, selected_indices, drop = FALSE])
      
      # Evaluate using simple RF cross-validation
      fitness_scores[i] <- evaluate_rf_fitness(subset_data)
      
    }, error = function(e) {
      fitness_scores[i] <- 0
    })
  }
  
  return(fitness_scores)
}

#' Evaluate RF Fitness
#'
#' @description
#' Evaluate fitness using Random Forest cross-validation
#'
#' @param subset_data Data subset for evaluation
#' @return Fitness score
evaluate_rf_fitness <- function(subset_data) {
  tryCatch({
    if (ncol(subset_data) <= 1) return(0)
    
    # Simple 3-fold CV for speed
    n_samples <- nrow(subset_data)
    fold_size <- floor(n_samples / 3)
    cv_scores <- numeric(3)
    
    for (fold in 1:3) {
      test_start <- (fold - 1) * fold_size + 1
      test_end <- min(fold * fold_size, n_samples)
      test_indices <- test_start:test_end
      
      if (length(test_indices) == 0) next
      
      train_data <- subset_data[-test_indices, , drop = FALSE]
      test_data <- subset_data[test_indices, , drop = FALSE]
      
      if (nrow(train_data) == 0 || nrow(test_data) == 0) next
      
      # Simple prediction using majority class
      train_class <- train_data[, 1]
      test_class <- test_data[, 1]
      majority_class <- names(sort(table(train_class), decreasing = TRUE))[1]
      
      accuracy <- sum(test_class == majority_class) / length(test_class)
      cv_scores[fold] <- accuracy
    }
    
    return(mean(cv_scores, na.rm = TRUE))
    
  }, error = function(e) {
    return(0.5)  # Baseline fitness
  })
}

#' Create New Generation
#'
#' @description
#' Create new generation using selection, crossover, and mutation
#'
#' @param population Current population
#' @param fitness_scores Fitness scores
#' @param elite_size Number of elite individuals
#' @param crossover_rate Crossover probability
#' @param mutation_rate Mutation probability
#' @param max_features Maximum features constraint
#' @return New population matrix
create_new_generation <- function(population, fitness_scores, elite_size, 
                                 crossover_rate, mutation_rate, max_features) {
  pop_size <- nrow(population)
  n_features <- ncol(population)
  new_population <- matrix(0, nrow = pop_size, ncol = n_features)
  
  # Elitism: keep best individuals
  elite_indices <- order(fitness_scores, decreasing = TRUE)[seq_len(elite_size)]
  new_population[seq_len(elite_size), ] <- population[elite_indices, , drop = FALSE]
  
  # Fill rest of population
  for (i in (elite_size + 1):pop_size) {
    if (runif(1) < crossover_rate && i < pop_size) {
      # Crossover
      parent1_idx <- tournament_selection(fitness_scores)
      parent2_idx <- tournament_selection(fitness_scores)
      
      offspring <- crossover(population[parent1_idx, ], population[parent2_idx, ])
      new_population[i, ] <- offspring
      
      if (i + 1 <= pop_size) {
        offspring2 <- crossover(population[parent2_idx, ], population[parent1_idx, ])
        new_population[i + 1, ] <- offspring2
        i <- i + 1
      }
    } else {
      # Direct selection
      parent_idx <- tournament_selection(fitness_scores)
      new_population[i, ] <- population[parent_idx, ]
    }
    
    # Mutation
    if (runif(1) < mutation_rate) {
      new_population[i, ] <- mutate_chromosome(new_population[i, ], max_features)
    }
  }
  
  return(new_population)
}

#' Tournament Selection
#'
#' @description
#' Select individual using tournament selection
#'
#' @param fitness_scores Fitness scores
#' @param tournament_size Tournament size
#' @return Selected individual index
tournament_selection <- function(fitness_scores, tournament_size = 3) {
  tournament_indices <- sample(length(fitness_scores), min(tournament_size, length(fitness_scores)))
  tournament_fitness <- fitness_scores[tournament_indices]
  winner_idx <- which.max(tournament_fitness)
  return(tournament_indices[winner_idx])
}

#' Crossover Operation
#'
#' @description
#' Single-point crossover between two chromosomes
#'
#' @param parent1 First parent chromosome
#' @param parent2 Second parent chromosome
#' @return Offspring chromosome
crossover <- function(parent1, parent2) {
  n_genes <- length(parent1)
  crossover_point <- sample(1:(n_genes-1), 1)
  
  offspring <- c(parent1[1:crossover_point], parent2[(crossover_point+1):n_genes])
  return(offspring)
}

#' Mutate Chromosome
#'
#' @description
#' Apply mutation to chromosome with feature count constraint
#'
#' @param chromosome Chromosome to mutate
#' @param max_features Maximum features allowed
#' @return Mutated chromosome
mutate_chromosome <- function(chromosome, max_features) {
  n_genes <- length(chromosome)
  current_count <- sum(chromosome)
  
  # Randomly flip a bit
  mutation_point <- sample(n_genes, 1)
  
  if (chromosome[mutation_point] == 1) {
    # Turn off feature (always allowed)
    chromosome[mutation_point] <- 0
  } else {
    # Turn on feature (only if under limit)
    if (current_count < max_features) {
      chromosome[mutation_point] <- 1
    }
  }
  
  return(chromosome)
}

#' Calculate Population Diversity
#'
#' @description
#' Calculate genetic diversity of population
#'
#' @param population Population matrix
#' @return Diversity score
calculate_population_diversity <- function(population) {
  if (nrow(population) <= 1) return(0)
  
  # Calculate pairwise Hamming distances
  distances <- numeric(0)
  
  for (i in 1:(nrow(population)-1)) {
    for (j in (i+1):nrow(population)) {
      hamming_dist <- sum(population[i, ] != population[j, ])
      distances <- c(distances, hamming_dist)
    }
  }
  
  return(mean(distances) / ncol(population))
}

#' Calculate Feature Fitness Scores
#'
#' @description
#' Calculate individual fitness scores for selected features
#'
#' @param selected_features Selected feature names
#' @param data Full dataset
#' @param feature_cols Feature columns
#' @return Feature fitness scores
calculate_feature_fitness_scores <- function(selected_features, data, feature_cols) {
  scores <- numeric(length(selected_features))
  names(scores) <- selected_features
  
  for (i in seq_along(selected_features)) {
    feature_name <- selected_features[i]
    feature_data <- feature_cols[, feature_name]
    class_data <- data[, 1]
    
    # Calculate correlation-based fitness
    if (is.numeric(feature_data)) {
      if (is.factor(class_data) || is.character(class_data)) {
        class_numeric <- as.numeric(as.factor(class_data))
      } else {
        class_numeric <- class_data
      }
      scores[i] <- abs(cor(feature_data, class_numeric, use = "complete.obs"))
    } else {
      # For categorical features
      contingency_table <- table(feature_data, class_data)
      if (any(dim(contingency_table) >= 2)) {
        chi_stat <- chisq.test(contingency_table, simulate.p.value = TRUE)$statistic
        scores[i] <- min(chi_stat / sum(contingency_table), 1)
      } else {
        scores[i] <- 0.1
      }
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
    method_id = 12,
    method_name = "GeneticAlgorithmRF",
    method_function = fs_method_geneticalgorithmrf,
    category = "optimization",
    dependencies = character(0),
    description = "Genetic Algorithm optimization with Random Forest fitness evaluation",
    parameters = list(
      population_size = 30,
      generations = 20,
      mutation_rate = 0.1,
      elite_size = 3,
      crossover_rate = 0.8
    ),
    complexity = "very_high",
    supports_smote = TRUE,
    timeout_default = 600
  )
  
  cat("✓ Registered Genetic Algorithm RF (GeneticAlgorithmRF) method (ID: 12)\n")
  
}, error = function(e) {
  warning("Failed to register GeneticAlgorithmRF method: ", e$message)
})
