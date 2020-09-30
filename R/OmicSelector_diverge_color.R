#' OmicSelector_diverge_color
#'
#' Helper for plotting.
#'
#' @export
OmicSelector_diverge_color <- function(data,centeredOn=0){
  suppressMessages(library(classInt))
  nHalf=50
  Min <- min(data,na.rm=TRUE)
  Max <- max(data,na.rm=TRUE)
  Thresh <- centeredOn
  pal<-colorRampPalette(c("blue", "white", "red"))(n = 11)
  rc1<-colorRampPalette(colors=c(pal[1],pal[2]),space="Lab")(10)
  for(i in 2:10){
    tmp<-colorRampPalette(colors=c(pal[i],pal[i+1]),space="Lab")(10)
    rc1<-c(rc1,tmp)
  }
  rb1 <- seq(Min, Thresh, length.out=nHalf+1)
  rb2 <- seq(Thresh, Max, length.out=nHalf+1)[-1]
  rampbreaks <- c(rb1, rb2)
  cuts <- classIntervals(data, style="fixed",fixedBreaks=rampbreaks)
  return(list(cuts,rc1))
}
