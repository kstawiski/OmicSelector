\name{os_make_benchmark_folds}
\alias{os_make_benchmark_folds}
\title{Reproduce the archived benchmark grouped folds}
\usage{os_make_benchmark_folds(y, group_id, k = 5L, seed = 42L)}
\arguments{
\item{y}{Numeric binary outcome vector (0 and 1).}
\item{group_id}{Nonmissing group identifiers, one per row.}
\item{k}{Number of folds, at least two.}
\item{seed}{Random seed.}}
\value{An unnamed list of sorted integer test-row indices.}
\description{Groups are enumerated in first-row order, classified by rounded mean outcome, shuffled within outcome (0 then 1), and assigned round-robin.}
\details{This reproduces the benchmark harness algorithm and changes the global RNG state through \code{set.seed(seed)}. A mixed-outcome group with mean 0.5 is assigned to class 0. Supply original row order and the same \code{RNGkind} settings (the archive uses the R default). The existing \code{\link{os_make_grouped_stratified_folds}} balances row counts and retains its different partitions and mixed-group checks.}
