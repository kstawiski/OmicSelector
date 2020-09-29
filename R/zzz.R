.onAttach <- function(libname, pkgname) {
  # options(rgl.useNULL = TRUE)
  # OmicSelector_setup(keras = FALSE, msg = FALSE)
  packageStartupMessage("\n\nWelcome to OmicSelector!\nAuthors: Konrad Stawiski M.D. (konrad@konsta.com.pl) and Marcin Kaszkowiak.\n\n\nYou can start with running OmicSelector_setup() to be sure that everything is properlly installed.\nFor more details go to https://kstawiski.github.io/OmicSelector/\n")
}

.onLoad <- function(libname, pkgname) {
  options(rgl.useNULL = TRUE)
}
