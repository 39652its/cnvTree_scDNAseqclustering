# 7.1: pqArmClustering() stands for pqArm clustering method
#' Perform pqArm clustering step in scDNA-seq data clustering Workflow
#'
#' This function executes the pqArm clustering step in the single-cell DNA sequencing (scDNA-seq)
#' data clustering workflow, grouping cells based on chromosomal arm-level copy number variations.
#'
#' @param input A named list where each element is a `GRanges` object representing a single cell.
#' @param pqArm_file In-build cytoband template for selection: `hg38`, `hg19`, `mm10`, `mm39`.
#' Or a filepath of a table for cytoband information seen on Giemsa-stained chromosomes.
#' It should include the following columns:
#'   - `chrom`: Reference sequence chromosome or scaffold.
#'   - `chromStart`: Start position in genoSeq.
#'   - `chromEnd`: End position in genoSeq.
#'   - `name`: Name of cytogenetic band.
#'   - `gieStain`: Giemsa stain results.
#'
#' @return A data frame recording the pqArm clustering results for each cell, tracking the clustering history at this stage.
#'
#' @export
#'
#' @examples
#' file_path <- system.file("extdata", "example_data.rds", package = "cnvTree")
#' Example_data <- changeFormat(file = file_path, core = 4)
#'
#' pqArm_result <- pqArmClustering(input = Example_data, pqArm_file = "hg38")
#'
pqArmClustering <- function(input, pqArm_file){
  message("=== Step 02: pqArm Clustering ===")

  Clustering_output <- data.frame(cluster = 1,
                                  cell = names(input))

  pqArm_result <- NULL
  for(Label in unique(Clustering_output$cluster)){
    ptm <- startTimed("pqArm Clustering ...")
    Smooth_pqCN <- pqArm_CN(input = input, Cluster_label = Label, Clustering_output = Clustering_output, pqArm_file = pqArm_file)

    pqArm_cluster <- pqArm_clustering(matrix = Smooth_pqCN, Label = Label)
    pqArm_cluster_summary <- pqArm_clustering_summary(matrix = pqArm_cluster, Label = Label)
    pqArm_cluster <- merge(pqArm_cluster, pqArm_cluster_summary)

    pqArm_result <- rbind(pqArm_result, pqArm_cluster)

    endTimed(ptm)
  }

  return(pqArm_result)
}


# 7.2: Reclustering() stands for Reclustering method
#' Perform Re-clustering in scDNA-seq data clustering workflow
#'
#' This function executes the re-clustering step in the single-cell DNA sequencing (scDNA-seq)
#' data clustering workflow, reclustering clusters based on differences in chromosomal arm-level
#' copy number variations.
#'
#' @param input A named list where each element is a `GRanges` object representing a single cell.
#' @param pqArm_output A data frame recording the pqArm clustering results for each cell.
#'   This table tracks the clustering history up to this stage.
#' @param pqArm_file In-build cytoband template for selection: `hg38`, `hg19`, `mm10`, `mm39`.
#' Or a filepath of a table for cytoband information seen on Giemsa-stained chromosomes.
#' It should include the following columns:
#'   - `chrom`: Reference sequence chromosome or scaffold.
#'   - `chromStart`: Start position in genoSeq.
#'   - `chromEnd`: End position in genoSeq.
#'   - `name`: Name of cytogenetic band.
#'   - `gieStain`: Giemsa stain results.
#' @param difratio_chr A numeric value defining the threshold for acceptable difference ratios across different chromosomes.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @return A data frame recording both pqArm clustering and re-clustering results for each cell,
#'   maintaining the clustering history across steps.
#'
#' @export
#'
#' @examples
#' file_path <- system.file("extdata", "example_data.rds", package = "cnvTree")
#' Example_data <- changeFormat(file = file_path, core = 4)
#'
#' pqArm_result <- pqArmClustering(input = Example_data, pqArm_file = "hg38")
#' Reclustering_output <- Reclustering(input = Example_data,
#'                                     pqArm_output = pqArm_result,
#'                                     pqArm_file = "hg38")
#'
Reclustering <- function(input, pqArm_output, pqArm_file, difratio_chr=0.3){
  message("=== Step 03: Reclustering ===")

  Recluster_Output <- NULL
  for(Label in unique(pqArm_output$cluster)){
    ptm <- startTimed("Reclustering clusters to make cluster bigger ...")
    pqArm_merge_target <- NULL
    New_pqArm_clustering <- pqArm_recluster(pqArm_cluster = pqArm_output, Cluster = Label)

    # potential merge target must more than one (with less10 clusters( >=2 & <10), or more than one more10 clusters(>2))
    if(length(New_pqArm_clustering)>1){
      pqArm_dif <- pqArm_reclustering_dif(input = input, pqArm_recluster_sim = New_pqArm_clustering,
                                          pqArm_cluster = pqArm_output, Cluster = Label, pqArm_file = pqArm_file)
      pqArm_merge_target <- pqArm_reclusterBy_ratio_target(pqArm_cluster = pqArm_output, Cluster = Label, pqReclsut_sim = pqArm_dif, difratio_chr = difratio_chr)
    }

    # is.null(pqArm_merge_target)
    if(is.null(pqArm_merge_target) == FALSE){
      Recluster_CellID <- pqArm_recluster_result(pqArm_cluster = pqArm_output, pqReclsut_target = pqArm_merge_target)
      Recluster_Output <- rbind(Recluster_Output, Recluster_CellID)
    } else {
      Recluster_CellID <- pqArm_output %>%
        dplyr::filter(.data$cluster %in% Label) %>%
        dplyr::mutate(Recluster_pattern = .data$pqArm_pattern,
                      Recluster_cellnum = .data$pqArm_cellnum,
                      Recluster_cluster = .data$pqArm_cluster)
      Recluster_Output <- rbind(Recluster_Output, Recluster_CellID)
    }
    endTimed(ptm)
  }

  return(Recluster_Output)
}


# 7.3: SubcloneClustering() stands for Subclone clsutering method
#' Perform subclone clustering in scDNA-seq data clustering workflow
#'
#' This function executes the subclone clustering step in the single-cell DNA sequencing (scDNA-seq) data clustering workflow,
#' refining clusters based on copy number variations (CNVs) within subpopulations of cells.
#'
#' @param input A named list where each element is a `GRanges` object representing a single cell.
#' @param Reclustering_output A data frame recording pqArm clustering and re-clustering results for each cell.
#'   This table tracks the clustering history up to this stage.
#' @param min_cell An integer specifying the minimum number of cells required for a cluster to be retained.
#' @param overlap_region An integer representing the genomic region where copy number frequently changes.
#' @param dif_ratio A numeric value defining the tolerance threshold for copy number differences
#'   between cells within a cluster.
#'
#' @return A list containing two data frames:
#'   - `final_cluster_output`: Records the clustering history from the pqArm, re-clustering, and subclone clustering steps.
#'   - `Subclone_CN`: Records each subclone's unique chromosome segment template and its copy number.
#'     It includes the following columns:
#'     - `chr`: Chromosome name (chr1, chr2, ...).
#'     - `start`: Start position of the segment.
#'     - `end`: End position of the segment.
#'     - `region`: The defined region index
#'     - `Subclone`: Subclone identifier.
#'     - `CN`: Copy number value of the segment.
#'
#' @export
#'
#' @examples
#' file_path <- system.file("extdata", "example_data.rds", package = "cnvTree")
#' Example_data <- changeFormat(file = file_path, core = 4)
#'
#' pqArm_result <- pqArmClustering(input = Example_data, pqArm_file = "hg38")
#' Reclustering_output <- Reclustering(input = Example_data,
#'                                     pqArm_output = pqArm_result,
#'                                     pqArm_file = "hg38")
#' Subclone_output <- SubcloneClustering(input = Example_data,
#'                                      Reclustering_output = Reclustering_output)
#'
SubcloneClustering <- function(input, Reclustering_output, min_cell=5, overlap_region=10**7, dif_ratio=0.2){
  message("=== Step 04: Subclone Clustering ===")

  Final_output <- list()
  Subclone_cluster <- NULL
  Subclone_CNr <- NULL
  cell_clustering <- NULL

  # from overlap_region to number of bins in coverage
  binsize <- as.data.frame(input[[1]]$bins)$width[1]
  overlap_bp <- round(overlap_region/binsize, digits = 0)

  # Select enough cell number Reclusters, avoiding keep filtering
  R <- Reclustering_output %>% dplyr::filter(.data$Recluster_cellnum >= min_cell)
  for (Recluster_label in unique(R$Recluster_cluster)){
    ptm <- startTimed("Subclone clustering for Recluster ", Recluster_label, " ... \n")

    breakpoints <- collect_cluster_bp(input = input, Clustering_output = R, Recluster_label = Recluster_label)
    Cell_num <- R %>% dplyr::filter(.data$Recluster_cluster == Recluster_label) %>% dplyr::pull(.data$Recluster_cellnum)
    bp <- output_bp_covers(Template = breakpoints, binsize = binsize, overlap=overlap_bp, overlap_times = Cell_num[1]*0.5) #overlap_times = 0.5*cell
    event_region <- bp_events(input = input, Template = bp, binsize = binsize)
    consensus_bp_template <- bp_region(event = event_region, binsize = binsize)
    cell_CNregion <- Region_CN(input = input, Reclustering_output = Reclustering_output,
                               Recluster_label = Recluster_label, events = consensus_bp_template)

    if(is.null(cell_clustering) == TRUE){
      cell_clustering <- Subclone_clustering(CN_incells_input= cell_CNregion, event_region= consensus_bp_template,
                                             dif_ratio = dif_ratio, Subclone_num = 0)
    } else {
      Subclone_num = max(Subclone_cluster$Subclone)
      cell_clustering <- Subclone_clustering(CN_incells_input= cell_CNregion, event_region= consensus_bp_template,
                                             dif_ratio = dif_ratio, Subclone_num = Subclone_num)
    }
    endTimed(ptm)

    Subclone_cluster <- rbind(Subclone_cluster, cell_clustering)

    # Output: 1. Copy number in each region in each subclone
    s_CN <- Subclone_CNregion(sep_region = consensus_bp_template, CN_region = cell_CNregion, each_subclone = cell_clustering,
                              min_cell = min_cell, output = "SubcloneRegionCN")
    Subclone_CNr <- rbind(Subclone_CNr, s_CN)

  }

  colnames(Subclone_cluster) <- c("Subclone_cluster", "cellID", "Subclone_cellnum")
  Final_output <- list(final_cluster_output = dplyr::left_join(Reclustering_output, Subclone_cluster, by = "cellID"),
                       Subclone_CN = Subclone_CNr)


  return(Final_output)
}


# 7.4: scDNA_Output() stands for outputting scDNA clustering final outputs
#' Generate final output in scDNA-seq data clustering workflow
#'
#' This function produces the final output for the single-cell DNA sequencing (scDNA-seq)
#' clustering workflow, integrating clustering results and generating visualizations.
#'
#' @param input A named list where each element is a `GRanges` object representing a single cell.
#' @param Summary A list containing two data frames:
#'   - `final_cluster_output`: Records the clustering history from the pqArm, re-clustering, and subclone clustering steps.
#'   - `Subclone_CN`: Records each subclone's unique chromosome segment template and its copy number.
#'     It includes the following columns:
#'     - `chr`: Chromosome name (chr1, chr2, ...).
#'     - `start`: Start position of the segment.
#'     - `end`: End position of the segment.
#'     - `region`: The defined region index
#'     - `Subclone`: Subclone identifier.
#'     - `CN`: Copy number value of the segment.
#' @param pqArm_file In-build cytoband template for selection: `hg38`, `hg19`, `mm10`, `mm39`.
#' Or a filepath of a table for cytoband information seen on Giemsa-stained chromosomes.
#' It should include the following columns:
#'   - `chrom`: Reference sequence chromosome or scaffold.
#'   - `chromStart`: Start position in genoSeq.
#'   - `chromEnd`: End position in genoSeq.
#'   - `name`: Name of cytogenetic band.
#'   - `gieStain`: Giemsa stain results.
#' @param output.dir A character string specifying the directory where the output files will be saved.
#' @param consecutive_region A numeric value defining the minimum length threshold for filtering copy number variation (CNV) regions.
#' @param cellcutoff A numeric value specifying the minimum number of cells required for a cluster to be retained.
#' @param sexchromosome_plot A logical value. If `TRUE`, the output includes plots with sex chromosome copy number information.
#' @param smoothheatmap A logical value. If `TRUE`, the output applies smoothing over a 10⁶ bp range in chromosome copy number visualizations.
#'
#' @return The function generates final clustering results and visualizations, saving them to the specified output directory.
#' @export
#'
#' @examples
#' \dontrun{
#' file_path <- system.file("extdata", "example_data.rds", package = "cnvTree")
#' Example_data <- changeFormat(file = file_path, core = 4)
#'
#' pqArm_result <- pqArmClustering(input = Example_data, pqArm_file = "hg38")
#' Reclustering_output <- Reclustering(input = Example_data,
#'                                     pqArm_output = pqArm_result,
#'                                     pqArm_file = "hg38")
#' Subclone_output <- SubcloneClustering(input = Example_data,
#'                                    Reclustering_output = Reclustering_output)
#' scDNA_Output(input = Example_data,
#'              Summary = Final_output,
#'              pqArm_file = "hg38")
#' }
#'
scDNA_Output <- function(input, Summary, pqArm_file, output.dir=getwd(), consecutive_region=10**7, cellcutoff=5, sexchromosome_plot=FALSE, smoothheatmap=TRUE){
  message("=== Step 05: Output cnvTree results ===")


  # 1. cellID summary
  writeOutput(data = Summary$final_cluster_output, filename = "/cnvTree.scDNAseq_grouping", path = output.dir)
  message("Output /cnvTree.scDNAseq_grouping.txt is done.")

  # 2. Each region copy number to subclone
  writeOutput(data = Summary$Subclone_CN, filename = "/cnvTree.scDNAseq_SubcloneRegionCN", path = output.dir)
  message("Output /cnvTree.scDNAseq_SubcloneRegionCN.txt is done.")

  # 3. All CNV region in the sample
  CNV_Data <- Total_cnvRegion(input = input, Template = Summary$Subclone_CN, pqArm_file = pqArm_file, consecutive_region = consecutive_region)
  CNV_Data <- cnvRegion.toPQarm(FILE = CNV_Data, pqArm_file = pqArm_file)
  writeOutput(data = CNV_Data, filename = "/cnvTree.scDNAseq_DefinedCNVregion", path = output.dir)
  message("Output /cnvTree.scDNAseq_DefinedCNVregion.txt is done.")

  # 4. DNA superimpose
  DNA_superimpose <- scDNA.superimpose(Template = Summary, DefinedCNVs = CNV_Data)
  DNA_superimpose <- scDNA.clustering(Template = DNA_superimpose)
  writeOutput(data = DNA_superimpose$DNA_cluster, filename = "/cnvTree.scDNAseq_DNAcluster", path = output.dir)
  message("Output /cnvTree.scDNAseq_DNAcluster.txt is done.")


  # 5. Final cluster heatmap
  Totalcluster_pdf(
    Input = input, Template = Summary$final_cluster_output,
    pqArm_file = pqArm_file, cellcutoff = cellcutoff, step = "Subclone", sexchromosome_plot = sexchromosome_plot,
    FILEpath = output.dir, FILEname = "/cnvTree.scDNAseq_Grouping_fig.pdf" )
  message("Output /cnvTree.scDNAseq_Grouping_fig.pdf is done.")

  # 6. cluster heatmap with dendrogram
  scDNA_CNVpattern(Input=input, final_cluster = Summary$final_cluster_output, cellcutoff = cellcutoff,
                   sexchromosome = sexchromosome_plot, smoothing = smoothheatmap,
                   pqArm_file=pqArm_file, FILEpath=output.dir, FILEname=paste0("/cnvTree.scDNAseq_heatmap.png"))
  message("Output /cnvTree.scDNAseq_heatmap.png is done.")
}


# 7.5 cnvTree_scDNAclustering() for whole scDNA cell clustering steps in 2 steps
#' One-Step scDNA-seq Cell Clustering Pipeline
#'
#' This function executes the entire single-cell DNA sequencing (scDNA-seq) clustering workflow in one step,
#' including pqArm clustering, re-clustering, subclone clustering, and final CNV-based output generation.
#'
#' @param input A named list where each element is a `GRanges` object representing a single cell.
#' @param pqArm_file In-build cytoband template for selection: `hg38`, `hg19`, `mm10`, `mm39`.
#' Or a filepath of a table for cytoband information seen on Giemsa-stained chromosomes.
#' It should include the following columns:
#'   - `chrom`: Reference sequence chromosome or scaffold.
#'   - `chromStart`: Start position in genoSeq.
#'   - `chromEnd`: End position in genoSeq.
#'   - `name`: Name of cytogenetic band.
#'   - `gieStain`: Giemsa stain results.
#' @param difratio_chr A numeric value defining the threshold for acceptable difference ratios across different chromosomes during re-clustering.
#' @param min_cell_subclone An integer specifying the minimum number of cells required for a cluster to be retained during subclone clustering.
#' @param overlap_region_subclone An integer representing genomic regions where copy number frequently changes in the subclone clustering step.
#' @param dif_ratio_subclone A numeric value defining the tolerance threshold for copy number differences between cells within a cluster in the subclone clustering step.
#' @param output.dir A character string specifying the directory where output files will be saved.
#' @param consecutive_region_output A numeric value defining the minimum CNV region length threshold for filtering CNV events in the final output.
#' @param cellcutoff_output A numeric value specifying the minimum number of cells required for a cluster to be retained in the final output.
#' @param sexchromosome_plot A logical value. If `TRUE`, the final output includes plots with sex chromosome copy number information.
#' @param smoothheatmap A logical value. If `TRUE`, the final output applies smoothing over a 10⁶ bp range in chromosome copy number visualizations.
#'
#' @return The function performs complete clustering and CNV analysis, saving final results and visualizations in the specified output directory.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' file_path <- system.file("extdata", "example_data.rds", package = "cnvTree")
#' Example_data <- changeFormat(file = file_path, core = 4)
#'
#' cnvTree_scDNAclustering(input=Example_data,
#'                         pqArm_file="hg38")
#' }
#'
cnvTree_scDNAclustering <- function(input, pqArm_file, difratio_chr=0.3,
                                    min_cell_subclone=5, overlap_region_subclone=10**7, dif_ratio_subclone=0.2,
                                    output.dir=getwd(), consecutive_region_output=10**7, cellcutoff_output=5, sexchromosome_plot=FALSE, smoothheatmap=TRUE){
  # step 01: pqArm clustering
  pqArm_result <- pqArmClustering(input = input, pqArm_file = pqArm_file)

  # step 02: Re-clustering
  Reclustering_output <- Reclustering(input = input, pqArm_output = pqArm_result, pqArm_file = pqArm_file, difratio_chr=difratio_chr)

  # step 03: Subclone clustering
  Subclone_output <- SubcloneClustering(input = input, Reclustering_output = Reclustering_output,
                                     min_cell=min_cell_subclone, overlap_region=overlap_region_subclone, dif_ratio=dif_ratio_subclone)

  # step 04: Output
  scDNA_Output(input = input,
               Summary = Subclone_output,
               output.dir = output.dir,
               pqArm_file = pqArm_file,
               consecutive_region = consecutive_region_output,
               cellcutoff = cellcutoff_output,
               sexchromosome_plot = sexchromosome_plot,
               smoothheatmap = smoothheatmap)
}








