.onAttach <- function(libname, pkgname) {
  # options(rgl.useNULL = TRUE)
  # OmicSelector_setup(keras = FALSE, msg = FALSE)
  packageStartupMessage("\n\nWelcome to OmicSelector!\nAuthors: Konrad Stawiski M.D. (konrad.stawiski@umed.lodz.pl) and Marcin Kaszkowiak.\n\nFor more details go to https://biostat.umed.pl/OmicSelector/\n")
}

.onLoad <- function(libname, pkgname) {
  options(rgl.useNULL = TRUE)
}
