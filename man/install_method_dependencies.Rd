\name{install_method_dependencies}
\alias{install_method_dependencies}
\title{Install dependencies for selected methods}
\usage{install_method_dependencies(methods = NULL, families = NULL, python = TRUE,
  dry_run = FALSE, envname = "omicselector", lib = .libPaths()[1L],
  repos = "https://cloud.r-project.org")}
\arguments{
\item{methods}{Optional method identifiers. NULL with NULL families selects all built-ins.}
\item{families}{Optional family letters; combined with methods as a union.}
\item{python}{Install Python modules and download required weights. FALSE installs only R requirements.}
\item{dry_run}{Return the plan without installation, Python imports or downloads.}
\item{envname}{Reticulate virtual environment name or path.}
\item{lib}{R library for newly installed packages.}
\item{repos}{CRAN repository.}
}
\value{A dependency plan, invisibly after installation.}
\description{Installs missing required R packages, prepares a reticulate virtual environment and obtains cached foundation checkpoints. Repeated installation retains existing R packages and Python packages satisfying the plan. Foundation module versions are pinned to the inspected environment. Other Python requirements use the selected package resolver. The installer does not request upgrades of requirements already satisfied.}
\details{When no Python interpreter is available, reticulate installs Miniconda before creating the environment. Use a fresh R session if another Python interpreter is initialized. Set RETICULATE_PYTHON to the returned environment's interpreter in later sessions. The selbal package comes from malucalle/selbal through remotes. Network access is needed for missing packages and weights. The installer does not alter model-fitting defaults or install optional alternative implementations.}
\examples{
install_method_dependencies("ai-tabpfn", dry_run = TRUE)
}
