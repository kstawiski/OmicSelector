  
options(warn = -1)
if(file.exists("task.log")) { file.remove("task.log") }
con <- file("task.log")
sink(con, append=TRUE)
sink(con, append=TRUE, type = "message")
library(OmicSelector)

try({
  current = 1
  max = parallel::detectCores()
  if(max < 3) { current = 0 }
  while(current > 0.5) { 
  load = strsplit(system("cat /proc/loadavg", intern = T)," ")
  current = as.numeric(load[[1]][1])/max
  if(current > 0.5) {
    cat(paste0("Current server load: ", round(current*100,2), "% exceeds the threshold of 50%. The job waiting for resources to start...\n")); Sys.sleep(15);
  } else { cat(paste0("Current server load: ", round(current*100,2), "%. The job is starting...\n")); }}
})

cat("Loading benchmark settings...\n")
mm = read.csv("selected_benchmark.csv")
m = as.character(mm$m) # which methods to check?
# mxnet = ifelse(readLines("var_mxnet.txt", warn = F) == "TRUE", TRUE, FALSE)
# if (length(mxnet) == 0) { mxnet = FALSE }
# search_iters_mxnet = as.numeric(readLines("var_search_iters_mxnet.txt", warn = F))
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
  cores = detectCores(),
  input_formulas = readRDS("featureselection_formulas_final.RDS"),
  output_file = "benchmark.csv",
  # mxnet = mxnet,
  # gpu = gpu,
  algorithms = m,
  holdout = holdout,
  stamp = "OmicSelector", OmicSelector_docker = T
)
cat("\n\nBenchmarking done. Moving to the analysis of best signature...\n\n")
rmarkdown::render("best_signature.Rmd")


cat("[OmicSelector: TASK COMPLETED]")
OmicSelector_log("[OmicSelector: TASK COMPLETED]","task.log")
sink() 
sink(type = "message")