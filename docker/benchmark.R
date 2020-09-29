  
options(warn = -1)
if(file.exists("task.log")) { file.remove("task.log") }
con <- file("task.log")
sink(con, append=TRUE)
sink(con, append=TRUE, type = "message")

library(OmicSelector)

cat("Loading benchmark settings...\n")
mm = read.csv("selected_benchmark.csv")
m = as.character(mm$m) # which methods to check?
mxnet = ifelse(readLines("var_mxnet.txt", warn = F) == "TRUE", TRUE, FALSE)
if (length(mxnet) == 0) { mxnet = FALSE }
search_iters_mxnet = as.numeric(readLines("var_search_iters_mxnet.txt", warn = F))
search_iters = as.numeric(readLines("var_search_iters.txt", warn = F))
holdout = ifelse(readLines("var_holdout.txt", warn = F) == "TRUE", TRUE, FALSE)
# gpu = tensorflow::tf$test$is_gpu_available()
gpu = F # force GPU to F

cat("Ok. Starting benchmark. This will take a while.. be patient. You can monitor this by checking CPU-load and temp/benchmark.csv for preliminary results.\n")
OmicSelector_benchmark(
  wd = getwd(),
  search_iters = search_iters,
  keras_epochs = 5000,
  keras_threads = floor(parallel::detectCores()/2),
  search_iters_mxnet = search_iters_mxnet,
  cores = detectCores() - 1,
  input_formulas = readRDS("featureselection_formulas_final.RDS"),
  output_file = "benchmark.csv",
  mxnet = mxnet,
  gpu = gpu,
  algorithms = m,
  holdout = holdout,
  stamp = "OmicSelector"
)
cat("\n\nBenchmarking done. Moving to the analysis of best signature...\n\n")
rmarkdown::render("best_signiture.Rmd")


cat("[OmicSelector: TASK COMPLETED]")
sink() 
sink(type = "message")