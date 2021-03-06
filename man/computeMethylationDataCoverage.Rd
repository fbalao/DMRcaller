% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/coverage.R
\name{computeMethylationDataCoverage}
\alias{computeMethylationDataCoverage}
\title{Compute methylation data coverage}
\usage{
computeMethylationDataCoverage(methylationData, regions = NULL,
  context = "CG", breaks = NULL, proportion = TRUE)
}
\arguments{
\item{methylationData}{the methylation data stored as a \code{\link{GRanges}}
object with four metadata columns (see \code{\link{methylationDataList}}).}

\item{regions}{a \code{\link{GRanges}} object with the regions where to
compute the coverage. If \code{NULL}, the coverage is computed genome-wide.}

\item{context}{the context in which the DMRs are computed (\code{"CG"},
\code{"CHG"} or \code{"CHH"}).}

\item{breaks}{a \code{numeric} vector specifing the different values for the
thresholds when computing the coverage.}

\item{proportion}{a \code{logical} value indicating whether to compute the
proportion (\code{TRUE}) or raw counts (\code{FALSE}).}
}
\value{
a \code{vector} with the proportion (or raw count) of cytosines that
have the number of reads higher or equal than the threshold values specified
in the \code{breaks} \code{vector}.
}
\description{
This function computes the coverage for bisulfite sequencing data. It
returns a \code{vector} with the proportion (or raw count) of cytosines that
have the number of reads higher or equal than a \code{vector} of specified
thresholds.
}
\examples{
# load the methylation data
data(methylationDataList)

# compute coverage in CG context
breaks <- c(1,5,10,15)
coverage_CG_wt <- computeMethylationDataCoverage(methylationDataList[["WT"]],
                 context="CG", breaks=breaks)
}
\author{
Nicolae Radu Zabet and Jonathan Michael Foonlan Tsang
}
\seealso{
\code{\link{plotMethylationDataCoverage}},
\code{\link{methylationDataList}}
}

