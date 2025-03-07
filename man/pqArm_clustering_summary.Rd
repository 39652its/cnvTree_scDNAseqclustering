% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scDNAclustering_03_pqArmClustering.R
\name{pqArm_clustering_summary}
\alias{pqArm_clustering_summary}
\title{Summarize pqArm clustering step results}
\usage{
pqArm_clustering_summary(matrix, Label)
}
\arguments{
\item{matrix}{A data frame recording the clustering results for each cell,
including the clustering history at each step.}

\item{Label}{An integer specifying the cluster to compute in the arm-level CNV analysis.}
}
\value{
A data frame summarizing the pqArm clustering step, containing the following columns:
\itemize{
\item \code{pqArm_pattern}: The identified copy number patterns at the arm level.
\item \code{pqArm_pattern_cellnum}: The number of cells associated with each pqArm pattern.
\item \code{cluster}: The assigned cluster for each pattern.
\item \code{pqArm_cluster}: The final cluster grouping based on pqArm patterns.
}
}
\description{
This function provides a summary of the pqArm clustering process,
detailing the distribution of arm-level copy number variation (CNV) patterns
across different clusters.
}
\keyword{internal}
