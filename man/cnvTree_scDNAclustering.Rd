% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scDNAclustering_07_TotalSteps.R
\name{cnvTree_scDNAclustering}
\alias{cnvTree_scDNAclustering}
\title{One-Step scDNA-seq Cell Clustering Pipeline}
\usage{
cnvTree_scDNAclustering(
  input,
  pqArm_file,
  difratio_chr = 0.3,
  min_cell_subclone = 5,
  overlap_region_subclone = 10^7,
  dif_ratio_subclone = 0.2,
  output.dir = getwd(),
  consecutive_region_output = 10^7,
  cellcutoff_output = 5,
  sexchromosome_plot = FALSE,
  smoothheatmap = TRUE
)
}
\arguments{
\item{input}{A named list where each element is a \code{GRanges} object representing a single cell.}

\item{pqArm_file}{In-build cytoband template for selection: \code{hg38}, \code{hg19}, \code{mm10}, \code{mm39}.
Or a filepath of a table for cytoband information seen on Giemsa-stained chromosomes.
It should include the following columns:
\itemize{
\item \code{chrom}: Reference sequence chromosome or scaffold.
\item \code{chromStart}: Start position in genoSeq.
\item \code{chromEnd}: End position in genoSeq.
\item \code{name}: Name of cytogenetic band.
\item \code{gieStain}: Giemsa stain results.
}}

\item{difratio_chr}{A numeric value defining the threshold for acceptable difference ratios across different chromosomes during re-clustering.}

\item{min_cell_subclone}{An integer specifying the minimum number of cells required for a cluster to be retained during subclone clustering.}

\item{overlap_region_subclone}{An integer representing genomic regions where copy number frequently changes in the subclone clustering step.}

\item{dif_ratio_subclone}{A numeric value defining the tolerance threshold for copy number differences between cells within a cluster in the subclone clustering step.}

\item{output.dir}{A character string specifying the directory where output files will be saved.}

\item{consecutive_region_output}{A numeric value defining the minimum CNV region length threshold for filtering CNV events in the final output.}

\item{cellcutoff_output}{A numeric value specifying the minimum number of cells required for a cluster to be retained in the final output.}

\item{sexchromosome_plot}{A logical value. If \code{TRUE}, the final output includes plots with sex chromosome copy number information.}

\item{smoothheatmap}{A logical value. If \code{TRUE}, the final output applies smoothing over a 10⁶ bp range in chromosome copy number visualizations.}
}
\value{
The function performs complete clustering and CNV analysis, saving final results and visualizations in the specified output directory.
}
\description{
This function executes the entire single-cell DNA sequencing (scDNA-seq) clustering workflow in one step,
including pqArm clustering, consolidating, subclustering, and final CNV-based output generation.
}
\examples{
\dontrun{
file_path <- system.file("extdata", "example_data.rds", package = "cnvTree")
Example_data <- changeFormat(file = file_path, core = 4)

cnvTree_scDNAclustering(input=Example_data,
                        pqArm_file="hg38")
}

}
