% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scDNAclustering_05_SubcloneClustering.R
\name{event_region.bin}
\alias{event_region.bin}
\title{Convert event region data to continuous bin-level sites}
\usage{
event_region.bin(input, Template)
}
\arguments{
\item{input}{A named list where each element is a \code{GRanges} object representing a single cell.}

\item{Template}{A data frame recording copy number variation (CNV) events for each chromosome.
This table will be updated with bin-level site information.}
}
\value{
A data frame \code{Template} with additional columns:
\itemize{
\item \code{bins_level_min}: The minimum bin-level position for each event region.
\item \code{bins_level_max}: The maximum bin-level position for each event region.
These columns enable further bin-level calculations.
}
}
\description{
This function transforms event region data from discrete chromosome sites
into continuous bin-level sites for further analysis.
}
\keyword{internal}
