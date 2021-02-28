  
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

# suppressMessages(library(foreach))
# suppressMessages(library(doParallel))
# suppressMessages(library(parallel))
# suppressMessages(library(doParallel))
# cl <- makePSOCKcluster(useXDR = TRUE, 3,  outfile="task.log") # We do not recommend using more than 5 threads, beacuse some of the methods inhereditly use multicore processing.
# registerDoParallel(cl)

mm = read.csv("selected_methods.csv")
m = as.numeric(mm$m) # which methods to check?

prefer_no_features = as.numeric(readLines("var_prefer_no_features.txt", warn = F))
max_iterations = as.numeric(readLines("var_max_iterations.txt", warn = F))
timeout_sec = as.numeric(readLines("var_timeout_sec.txt", warn = F))

for (i in m) {
try(OmicSelector_OmicSelector(m = i, max_iterations = max_iterations, stamp = "fs", debug = F, # we set debug to false (to make the package smaller), you may also want to change stamp to something meaningful, max_iterations was set to 1 to recude the computational time.. in real life scenarios it is resonable to use at least 10 iterations.
                  prefer_no_features = prefer_no_features, # Few methods are filter rather than wrapper methods, thus requires the maximum number of maximum features.    
                  timeout_sec = timeout_sec))
}


# foreach(i = m) %do%
#for(i in m)
# {
#   suppressMessages(library(OmicSelector))
#   prefer_no_features = readLines("var_prefer_no_features.txt", warn = F)
#   max_iterations = readLines("var_max_iterations.txt", warn = F)
#   timeout_sec = readLines("var_timeout_sec.txt", warn = F)
#   # setwd("/OmicSelector/OmicSelector/vignettes") # change it you to your working directory
#   OmicSelector_OmicSelector(m = i, max_iterations = max_iterations, stamp = "fs", debug = F, # we set debug to false (to make the package smaller), you may also want to change stamp to something meaningful, max_iterations was set to 1 to recude the computational time.. in real life scenarios it is resonable to use at least 10 iterations.
#                   prefer_no_features = prefer_no_features, # Few methods are filter rather than wrapper methods, thus requires the maximum number of maximum features.    
#                   timeout_sec = timeout_sec) # We don't want to wait eternity in this tutorial, just 10 minutes. Timeout is useful for complicated methods. Depending on your CPU 2 days may be reasonable for larger projects. Note that some methods cannot be controled with timeout parameter.
# }
# stopCluster(cl)

selected_sets_of_miRNAs = OmicSelector_merge_formulas(max_miRNAs = prefer_no_features)

cat("[OmicSelector: TASK COMPLETED]")
OmicSelector_log("[OmicSelector: TASK COMPLETED]","task.log")
sink() 
sink(type = "message")
cat("[OmicSelector: TASK COMPLETED]")