options(repos=c(CRAN='https://cloud.r-project.org'),Ncpus=4L)
install.packages('https://cran.r-project.org/src/contrib/Archive/remotes/remotes_2.5.0.tar.gz',repos=NULL,type='source')
lock<-read.delim('/opt/locks/r-packages.lock.tsv',stringsAsFactors=FALSE)
# Resolve dependencies, then replace every locked package with its exact version.
for(i in seq_len(nrow(lock)))remotes::install_version(lock$Package[i],version=lock$Version[i],upgrade='never',dependencies=NA)
for(i in seq_len(nrow(lock)))stopifnot(as.character(packageVersion(lock$Package[i]))==lock$Version[i])
