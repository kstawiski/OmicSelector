suppressMessages(library(OmicSelector))
OmicSelector_benchmark(search_iters = 5, # 5 random hyperparameter sets will be checked; 5 is set here for speed purposes.. for real projects use more, like 5000...
            algorithms = c("mlp", "mlpML", "svmRadial", "svmLinear", "rf", "C5.0", "rpart",
                            "rpart2", "ctree"), # default set of methods, note that logistic regression (glm) is always included
            output_file = paste0("benchmark.csv")) # the main output
