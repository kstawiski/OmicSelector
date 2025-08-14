.onAttach <- function(libname, pkgname) {
  packageStartupMessage("\n",
    "\u2728 Welcome to OmicSelector v", utils::packageVersion("OmicSelector"), "!\n",
    "\ud83d\udcca Biomarker Signature Selection & Deep Learning\n",
    "\ud83d\ude80 Author: Konrad Stawiski M.D., Ph.D.\n",
    "\ud83d\udd17 More info: https://biostat.umed.pl/OmicSelector/\n",
    "\ud83d\udd0d Get started with: ?omics_select\n"
  )
}

.onLoad <- function(libname, pkgname) {
  # Set options for better compatibility
  options(rgl.useNULL = TRUE)
  
  # Suppress rgl graphics on headless systems
  Sys.setenv("RGL_USE_NULL" = "TRUE")
  
  # Initialize logging if logger package is available
  if (requireNamespace("logger", quietly = TRUE)) {
    logger::log_threshold("INFO", namespace = "OmicSelector")
  }
}

.onUnload <- function(libpath) {
  # Clean up any parallel backends
  if (requireNamespace("future", quietly = TRUE)) {
    future::plan(future::sequential)
  }
  
  # Clean up global cluster if exists
  if (exists(".omics_cluster", envir = .GlobalEnv)) {
    cluster <- get(".omics_cluster", envir = .GlobalEnv)
    try(parallel::stopCluster(cluster), silent = TRUE)
    rm(".omics_cluster", envir = .GlobalEnv)
  }
}
