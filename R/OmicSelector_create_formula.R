#' OmicSelector_create_formula
#'
#' Helper function to create formula based on selected miRNAs.
#'
#' @param selected_features Selected miRNAs as characted vector.
#' @return Formula "Class ~ ..." to be used in another functions.
#'
#' @export
OmicSelector_create_formula = function(selected_features) {
  selected_features<-selected_features[!is.na(selected_features)]
  as.formula(paste0("Class ~ ",paste0(as.character(selected_features), collapse = " + ")))
}
