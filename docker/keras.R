#!/usr/bin/env Rscript
r = getOption("repos")
r["CRAN"] = "https://cran.r-project.org"
options(repos = r)

# install.packages('devtools')
# devtools::install_github('rstudio/keras', upgrade = 'never')
install.packages(c('reticulate','opencv'))
reticulate::use_python('/opt/conda/bin/python')

require(tensorflow)
require(reticulate)
require(keras)

is_keras_available()
system('which python')
Sys.setenv(TENSORFLOW_PYTHON='/opt/conda/bin/python')
use_python('/opt/conda/bin/python')

py_discover_config('tensorflow')
py_discover_config('keras')
is_keras_available()

packages = c("kerasR",
             "ggplot2", 
             "dplyr",
             "magrittr",
             "zeallot",
             "tfruns")

for (package in packages) {
    install.packages(package)
}
