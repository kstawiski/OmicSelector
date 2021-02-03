## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, message=FALSE, warning=FALSE,
  comment = "#>"
)
knitr::opts_chunk$set(fig.width=12, fig.height=8)
knitr::opts_chunk$set(tidy.opts=list(width.cutoff=150),tidy=TRUE)
options(rgl.useNULL = TRUE)
options(warn=-1)
suppressMessages(library(dplyr))
set.seed(1)
options(knitr.table.format = "html")
library(OmicSelector)

