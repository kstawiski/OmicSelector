r = getOption("repos")
r["CRAN"] = "https://cran.r-project.org"
options(repos = r)

install.packages(c('repr', 'IRdisplay', 'evaluate', 'crayon', 'pbdZMQ', 'devtools', 'uuid', 'digest'))
devtools::install_github("IRkernel/IRkernel")
IRkernel::installspec(displayname = 'OmicSelector')