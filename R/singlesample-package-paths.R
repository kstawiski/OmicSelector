# Internal path helpers for promoted Paper 3 runner utilities.

.singlesample_extdata_path <- function(filename, supplied = NULL, what = filename,
                                 required = TRUE) {
  if (!is.null(supplied)) {
    supplied <- path.expand(as.character(supplied)[1L])
    if (!file.exists(supplied)) {
      stop(what, " not found at supplied path: ", supplied, call. = FALSE)
    }
    return(supplied)
  }

  ext <- system.file("extdata", filename, package = "OmicSelector")
  if (nzchar(ext) && file.exists(ext)) {
    return(ext)
  }

  if (isTRUE(required)) {
    stop(
      what, " is not shipped with this package build. Supply an explicit path ",
      "to the corresponding Paper 3 input.",
      call. = FALSE
    )
  }
  NULL
}

.singlesample_default_output <- function(output_dir, filename) {
  file.path(path.expand(as.character(output_dir)[1L]), filename)
}
