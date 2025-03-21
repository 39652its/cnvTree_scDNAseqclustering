% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scDNAclustering_03_pqArmClustering.R
\name{pqArm_file.cen}
\alias{pqArm_file.cen}
\title{Process UCSC cytoband data for "acen" and "gvar" regions}
\usage{
pqArm_file.cen(FILE)
}
\arguments{
\item{FILE}{A character string specifying the file path to the UCSC cytoband data file.}
}
\value{
A data frame containing labeled genomic ranges for "acen" and "gvar" cytoband
regions, which can be used for masking CNVs in experimental analyses.
}
\description{
This function extracts and formats cytoband information from the UCSC database,
specifically for the "acen" (centromeric) and "gvar" (variable heterochromatic)
cytoband types. The output is structured for defining masking ranges
to exclude copy number variations (CNVs) from downstream analyses.
}
\keyword{internal}
