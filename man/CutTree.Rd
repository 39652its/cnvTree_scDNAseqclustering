% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scDNAclustering_02_Clustering.R
\name{CutTree}
\alias{CutTree}
\title{Divide cells into two groups based on copy number variation in a specific cluster}
\usage{
CutTree(input, Template, Cluster_label)
}
\arguments{
\item{input}{A named list where each element is a \code{GRanges} object representing a single cell.}

\item{Template}{A data frame containing two columns:
\itemize{
\item \code{cellID}: Unique identifier for each cell.
\item \code{cluster}: Cluster assignment for each cell.
}}

\item{Cluster_label}{An integer specifying the cluster label used to divide the cells into two groups.}
}
\value{
A data frame with two columns:
\itemize{
\item \code{cellID}: The unique identifier of each cell.
\item \code{cluster}: The updated cluster assignment (k = 2).
}
}
\description{
This function separates cells into two groups based on copy number variation in a specified cluster label
from a provided clustering result table.
}
\keyword{internal}
