\name{method_requirements}
\alias{method_requirements}
\title{Requirements and measured costs for registered methods}
\usage{method_requirements(methods = NULL, families = NULL)}
\arguments{
\item{methods}{Optional method identifiers. NULL selects all unless families is set.}
\item{families}{Optional roster family letters. Methods and families select their union.}
}
\value{A data frame with one row per method. Dependencies are separated by semicolons. Optional R implementations are listed separately. Cost columns are medians of completed CPU cells. Missing measurements are NA. The \code{measurements} attribute contains full cells with source and hardware. User callbacks without dependency declarations have \code{requirements_known = FALSE}.}
\description{Lists dependencies beyond core Imports and, when available, measured PB1 fit time, per-specimen scoring time, peak memory and frozen model size. Measurements describe the recorded data, hardware and settings, not a speed guarantee.}
\examples{
method_requirements("ws-rclr-trimmed")
}
