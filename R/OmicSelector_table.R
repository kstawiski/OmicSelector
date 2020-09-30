#' OmicSelector_table
#'
#' A wrapper to draw pretty tables in rmarkdown document
#'
#' @param table Table to draw.
#' @param hight High (default: 400px)
#'
#' @export
OmicSelector_table = function(table, height = "400px", ...)
{
  suppressMessages(library(knitr))
  suppressMessages(library(rmarkdown))
  suppressMessages(library(kableExtra))
  if (is.null(sessionInfo()$loadedOnly$IRdisplay)) { # czy jestem w Jupyterze?
  if(nrow(table) >= 6 && ncol(table) >= 6) {
kable(table, "html", ...) %>%
    kable_styling() %>%
    scroll_box(width = "100%", height = height)
  } else { kable(table, ...) }
     } else {
      as.data.frame(table)
    }
}
