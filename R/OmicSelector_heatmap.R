#' OmicSelector_heatmap
#'
#' Draw a heatmap of selected miRNAs.
#'
#' @param x Matrix of log-transformed TPM-normalized counts with miRNAs in columns and cases in rows.
#' @param rlab Data frame of factors to be marked on heatmap (like batch or class). Maximum of 2 levels for every variable is supported.
#' @param zscore Whether to z-score values before clustering and plotting.
#' @param expression_name What should be written on the plot?
#' @param trim_min Trim lower than.. Useful for setting appropriate scale
#' @param trim_max Trim greater than.. Useful for setting appropriate scale
#' @param centered_on On which value should the scale be centered? If null - median will be used.
#' @param legend_pos Where should the legend should be? Default: topright
#' @param legend_cex How large should the legend should be? Default: 0.8
#'
#' @return Heatmap.
#'
#' @export
OmicSelector_heatmap = function(x = trainx[,1:10], rlab = data.frame(Batch = dane$Batch, Class = dane$Class), zscore = F, margins = c(10,10), expression_name = "log10(TPM)", trim_min = NULL, trim_max = NULL, centered_on = NULL, legend_pos = "topright", legend_cex = 0.8, ...) {
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
  assigcolor = character()
  assigcode = character()
  kolory = rep(palette(),20)[-1]
  kolor_i = 1
  for (i in 1:ncol(rlab)){
    rlab[,i] = as.factor(rlab[,i])
    o_ile = as.numeric(length(unique(rlab[,i])))
    #assigcode = c(assigcode, as.character(unique(rlab[,i]))) -> bug
    assigcode = c(assigcode, as.character(rev(levels(rlab[, i]))))
    assigcolor = c(assigcolor, kolory[kolor_i:(kolor_i+o_ile-1)])
    #levels(rlab[,i]) = topo.colors(length(unique(rlab[,i])))
    levels(rlab[,i]) = kolory[kolor_i:(kolor_i+o_ile-1)]
    kolor_i = kolor_i+o_ile
  }
  assig = data.frame(assigcode, assigcolor)
  #levels(rlab$Batch) = rainbow(length(unique(dane$Batch)))
  #levels(rlab$Class) = c("red","green") # red - cancer, green - control
  x2 = as.matrix(x)
  colnames(x2) = gsub("\\.","-", colnames(x2))



  if(zscore == F) {
    if(!is.null(trim_min)) { x2[x2<trim_min] = trim_min }
    if(!is.null(trim_max)) { x2[x2>trim_max] = trim_max }


    if(!is.null(centered_on)) { brks<-OmicSelector_diverge_color(x2, centeredOn = centered_on) }
    else { brks<-OmicSelector_diverge_color(x2, centeredOn = median(x2)) }



    # colors = seq(min(x2), max(x2), by = 0.01)
    # my_palette <- colorRampPalette(c("blue", "white", "red"))(n = length(colors) - 1)

    rlab = as.matrix(rlab)
    OmicSelector_heatmap.3(x2, hclustfun=OmicSelector_myclust, distfun=OmicSelector_mydist,
                 RowSideColors=t(rlab),
                 margins = margins,
                 KeyValueName=expression_name,
                 symm=F,symkey=F,symbreaks=T, scale="none",
                 col=as.character(brks[[2]]),
                 breaks=as.numeric(brks[[1]]$brks), assigcode=assigcode, assigcolor=assigcolor, legend_pos = legend_pos, legend_cex = legend_cex, ...
                 #legend = T
                 #,scale="column"
    )
  } else {
    x3 = x2
    for(i in 1:ncol(x2)) {
      x3[,i] = scale(x2[,i])
    }


    if(!is.null(trim_min)) { x3[x3<trim_min] = trim_min }
    if(!is.null(trim_max)) { x3[x3>trim_max] = trim_max }

    if(!is.null(centered_on)) { brks<-OmicSelector_diverge_color(x3, centeredOn = centered_on) }
    else { brks<-OmicSelector_diverge_color(x3, centeredOn = median(x3)) }

    # colors = seq(min(x2), max(x2), by = 0.01)
    # my_palette <- colorRampPalette(c("blue", "white", "red"))(n = length(colors) - 1)

    rlab = as.matrix(rlab)
    OmicSelector_heatmap.3(x3, hclustfun=OmicSelector_myclust, distfun=OmicSelector_mydist,
                 RowSideColors=t(rlab),
                 margins = margins,
                 KeyValueName=paste0("Z-score ",expression_name),
                 symm=F,symkey=F,symbreaks=T, scale="none",
                 col=as.character(brks[[2]]),
                 breaks=as.numeric(brks[[1]]$brks), assigcode=assigcode, assigcolor=assigcolor, legend_pos = legend_pos, legend_cex = legend_cex, ...
                 #legend = T
                 #,scale="column"
    )
  }
}
