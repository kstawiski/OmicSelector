#' OmicSelector_myclust
#'
#' Helper in heatmap creation. Which method of cultering should be used?
#'
#' @export
OmicSelector_myclust=function(c) {hclust(c,method="average")}
