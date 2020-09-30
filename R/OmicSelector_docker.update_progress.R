#' OmicSelector_docker.update_progress
#'
#' Helper for docker. Updates progress of analysis by 1 step.
#'
#' @export
#'
OmicSelector_docker.update_progress = function(verbose = F) {
    if (file.exists("var_progress.txt")) {
       old = readLines("var_progress.txt")
        new = as.numeric(old) + 1
        writeLines(as.character(new), "var_progress.txt", sep = "")
    } else {
        writeLines("1", "var_progress.txt")
    }
    if(verbose) { cat(paste0("Current progress: ", readLines("var_progress.txt") )) }
}
