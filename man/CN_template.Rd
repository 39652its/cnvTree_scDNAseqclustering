% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scDNAclustering_03_pqArmClustering.R
\name{CN_template}
\alias{CN_template}
\title{Generate a copy number segment template based on chromosomal arms}
\usage{
CN_template(input, pqArm_file)
}
\arguments{
\item{input}{A named list where each element is a \code{GRanges} object representing a single cell.}

\item{pqArm_file}{A table for cytoband information seen on Giemsa-stained chromosomes.
It should include the following columns:
\itemize{
\item \verb{chrom"}: Reference sequence chromosome or scaffold.
\item \code{chromStart}: Start position in genoSeq.
\item \code{chromEnd}: End position in genoSeq.
\item \code{name}: Name of cytogenetic band.
\item \code{gieStain}: Giemsa stain results.
}}
}
\value{
A data frame containing the defined genomic ranges for p/q arms across all chromosomes.
}
\description{
This function constructs a template for copy number segmentation
using cytoband information from Giemsa-stained chromosomes.
}
\keyword{internal}
