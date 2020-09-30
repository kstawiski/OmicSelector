#' OmicSelector_load_extension
#' 
#' Function for loading extension.
#' 
#' @param name Name of extension to be loaded.
#' 
#' @export 
OmicSelector_load_extension = function(name = "deeplearning") {
    if(file.exists(paste0(name, ".R"))) {
        source(paste0(name, ".R"))
    } else {
        link = paste0("https://raw.githubusercontent.com/kstawiski/OmicSelector/master/extensions/", name, ".R")
        devtools::source_url(link) 
        } 
}