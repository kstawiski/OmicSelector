.onAttach <- function(libname, pkgname) {
  # options(rgl.useNULL = TRUE)
  # OmicSelector_setup(keras = FALSE, msg = FALSE)
  packageStartupMessage("\n\nWelcome to OmicSelector!\nAuthor: Konrad Stawiski M.D., Ph.D. (konrad.stawiski@umed.lodz.pl)\n\nFor more details go to https://biostat.umed.pl/OmicSelector/\n")
}

.onLoad <- function(libname, pkgname) {
  options(rgl.useNULL = TRUE)
}
