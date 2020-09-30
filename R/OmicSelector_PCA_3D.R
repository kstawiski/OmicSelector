#' OmicSelector_PCA_3D
#'
#' Conduct PCA and create 3D scatterplot form first 3 components.
#'
#' @param ttpm_features Normalizaed counts in `ttpm` format.
#' @param meta Factor of cases labels that should be visualized on biplot.
#'
#' @return Plotly 3D object.
#'
#' @export
OmicSelector_PCA_3D = function(ttpm_features, meta) {
  suppressMessages(library(plyr))
  suppressMessages(library(dplyr))
  suppressMessages(library(edgeR))
  suppressMessages(library(epiDisplay))
  suppressMessages(library(rsq))
  suppressMessages(library(MASS))
  suppressMessages(library(Biocomb))
  suppressMessages(library(caret))
  suppressMessages(library(dplyr))
  suppressMessages(library(epiDisplay))
  suppressMessages(library(pROC))
  suppressMessages(library(ggplot2))
  suppressMessages(library(DMwR))
  suppressMessages(library(ROSE))
  suppressMessages(library(gridExtra))
  suppressMessages(library(gplots))
  suppressMessages(library(devtools))
  suppressMessages(library(stringr))
  suppressMessages(library(data.table))
  suppressMessages(library(tidyverse))
  dane.pca <- prcomp(ttpm_features, scale. = TRUE)
  suppressMessages(library(plotly))

  if(!is.null(sessionInfo()$loadedOnly$IRdisplay)) {
    message("Running this code in the Jupyter enviorment may cause trouble like crash of R kernel.")
    invisible(readline(prompt="Press [Enter] to continue or [Ctrl+C] to cancel."))
  }

  for(i in colnames(ttpm_features)) {
    if(!is.numeric(ttpm_features[, i])) {
      stop("Please provide a dataframe with only numeric variables")
    }
  }

  if(is.data.frame(meta)) {
    stop("Please provide a single categorical vector")
  }

  pc = as.data.frame(dane.pca$x)
  pc = cbind(pc, meta)

  p <- plot_ly(data = pc, x = ~PC1, y = ~PC2, z = ~PC3, color = ~meta) %>%
    add_markers() %>%
    layout(scene = list(xaxis = list(title = 'PC1'),
                        yaxis = list(title = 'PC2'),
                        zaxis = list(title = 'PC3')))

  p
}
