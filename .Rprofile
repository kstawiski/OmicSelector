# OmicSelector 2.0 R Profile
# Activates renv for reproducible package management

# Activate renv if available
if (file.exists("renv/activate.R")) {
  source("renv/activate.R")
}

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Prevent accidental global library usage
options(renv.config.auto.snapshot = FALSE)

# Informative startup message
if (interactive()) {
  message("OmicSelector 2.0 Development Environment")
  message("Using renv for package management. Run `renv::status()` to check.")
}
