\name{available_methods}
\alias{available_methods}
\title{Inspect installed method dependencies}
\usage{available_methods(methods = NULL, families = NULL, python = TRUE, weights = TRUE)}
\arguments{
\item{methods}{Optional method identifiers.}
\item{families}{Optional roster family letters; combined with methods as a union.}
\item{python}{Check imports in the active reticulate environment. FALSE leaves Python methods unverified.}
\item{weights}{Check local foundation weights without downloading. FALSE leaves weight requirements unverified.}
}
\value{Requirements plus available, missing_r, missing_python, missing_weights, reason and install_command. Unchecked or undeclared requirements give NA availability.}
\description{Checks installed R packages, importable Python modules and cached weights. Does not install or download. Imports initialize the configured Python interpreter. Dependency availability does not establish that a particular panel or metadata set satisfies the existing method audits.}
\examples{
available_methods("ws-rclr-trimmed", python = FALSE)
}
