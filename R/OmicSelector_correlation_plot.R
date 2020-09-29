#' OmicSelector_correlation_plot
#'
#' Draw a correlation plot with reference f(x)=x axis and correlation metrics.
#'
#' @param var1 First numeric variable.
#' @param var2 Second numeric variable.
#' @param labvar1 Label to be written on plot for fist variable.
#' @param labvar2 Label to be written on plot for second variable.
#' @param yx Logical varable. TRUE if you want f(x)=x axis for reference (usefull for calibration plots).
#' @param gdzie_legenda Where the legend should be placed?
#'
#' @export
OmicSelector_correlation_plot = function(var1, var2, labvar1, labvar2, title, yx = T, metoda = 'pearson', gdzie_legenda = "topleft") {
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
  var1 = as.numeric(as.character(var1))
  var2 = as.numeric(as.character(var2))
  plot(var1, var2, pch = 19, cex=0.5, xlab=labvar1, ylab=labvar2, main=title)
  if (yx == T) { abline(0,1, col='gray') }
  abline(fit <- lm(var2 ~ var1), col='black', lty = 2)
  temp = cor.test(var1, var2, method = metoda)
  if (metoda=='pearson') {
    if (temp$p.value < 0.0001) {
      legend(gdzie_legenda, bty="n", legend=paste("r =",round(cor(var1,var2,use="complete.obs"),2),", p < 0.0001\nadj. R2 =", format(summary(fit)$adj.r.squared, digits=4)))
    } else {
      legend(gdzie_legenda, bty="n", legend=paste("r =",round(cor(var1,var2,use="complete.obs"),2),", p =",round(temp$p.value, 4),"\nadj. R2 =", format(summary(fit)$adj.r.squared, digits=4))) } }
  if (metoda=="spearman")
  { if (temp$p.value < 0.0001) {
    legend(gdzie_legenda, bty="n", legend=paste("rho =",round(cor(var1,var2,use="complete.obs"),2),", p < 0.0001"))
  } else {
    legend(gdzie_legenda, bty="n", legend=paste("rho =",round(cor(var1,var2,use="complete.obs"),2),", p =",round(temp$p.value, 4))) } }

}
