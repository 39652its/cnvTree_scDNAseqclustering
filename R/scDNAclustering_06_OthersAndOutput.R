# 6.1: startTimed(), endTimed() for time calculating
#' Record current time for function started timing
#'
#' This utility function records the current system time, typically used
#' for timing the execution of other functions. It is useful for benchmarking
#' or logging the duration of function calls.
#'
#' @param ... Additional arguments passed to methods.
#' Currently not used but included for compatibility and extensibility.
#'
#' @return A POSIXct object representing the current system time at the moment this function is called.
#'
#' @keywords internal
#'
startTimed <- function(...){
  x <- paste0(..., collapse = "")
  message(x, appendLF = FALSE)
  ptm <- proc.time()
  return(ptm)
}


#' Calculate time elapsed since start time
#'
#' This utility function calculates the time elapsed since a recorded start time,
#' typically used for measuring function execution duration. It provides a message
#' displaying the time consumed between two lines of code.
#'
#' @param ptm A POSIXct object representing the start time, typically obtained from a call to \code{record_time()}.
#'
#' @return A message showing the time consumed between the start time and the moment this function is called.
#'
#' @keywords internal
#'
endTimed <- function(ptm){
  time <- proc.time() - ptm
  message(" ", round(time[3], 2), "s")
}

# 6.3: GenomeHeatmap() function for plotting CN pattern in each cells
#' Plot copy number pattern for each cell in heatmap
#'
#' This function visualizes the copy number variation (CNV) pattern for each cell,
#' using cytoband information to annotate chromosomal regions. It supports optional
#' inclusion of sex chromosome CNV data in the output plot.
#'
#' @param Input A named list where each element is a `GRanges` object representing a single cell.
#' @param cellID A character vector specifying the `cellID`s of the cells for heatmap plotting.
#' @param pqArm_file A table for cytoband information seen on Giemsa-stained chromosomes.
#' It should include the following columns:
#'   - `chrom`: Reference sequence chromosome or scaffold.
#'   - `chromStart`: Start position in genoSeq.
#'   - `chromEnd`: End position in genoSeq.
#'   - `name`: Name of cytogenetic band.
#'   - `gieStain`: Giemsa stain results.
#' @param sexchromosome Logical. If `TRUE`, the plot includes copy number information for sex chromosomes. Defaults to `FALSE`.
#'
#' @return A `ggplot` object displaying the copy number pattern for each cell.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' file_path <- system.file("extdata", "example_data.rds", package = "cnvTree")
#' Example_data <- changeFormat(file = file_path, core = 4)
#'
#' GenomeHeatmap(
#'   Input = Example_data,
#'   cellID = names(Example_data)[1:10],
#'   pqArm_file = "hg38",
#'   sexchromosome = TRUE
#' )
#'
#'
GenomeHeatmap <- function(Input, cellID, pqArm_file, sexchromosome=FALSE){
  # import CN template and total cell copy number matrix
  CN_bins_template <- CN_template(input = Input, pqArm_file = pqArm_file)
  CN_chr_template <- CN_bins_template %>%
    dplyr::group_by(.data$chr) %>%
    dplyr::summarise(length = max(.data$end) - min(.data$start)) %>%
    dplyr::mutate(X_cum = c(0, cumsum(as.numeric(length[-length(length)]))))

  # Import the HC clustering cell order
  cellOrder <- clusterbyHMM(input = Input, selected = cellID)
  cellOrder <- cellOrder$IDorder %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "cellID")
  Input <- Input[cellOrder$cellID]

  # Import copy number data
  numofcell <- length(cellID)
  mat_input <- lapply(1:numofcell, function(i) segment_transform(data = Input[[i]], index = i, CN_chr_template = CN_chr_template))
  mat_input <- dplyr::bind_rows(mat_input)

  # chrX, chrY in plotting output is optional
  if(sexchromosome==FALSE){
    CN_bins_template <- CN_bins_template %>% dplyr::filter(!.data$chr %in% c("chrX", "chrY"))
    CN_chr_template <- CN_chr_template %>% dplyr::filter(!.data$chr %in% c("chrX", "chrY"))
    mat_input <- mat_input %>% dplyr::filter(!.data$chr %in% c("chrX", "chrY"))
  }

  # mat_input$CN[which(mat_input$CN>5)] <- 5

  # plot chr site
  Chr_tmp <- CN_chr_template %>%
    dplyr::mutate(chr = sub("^chr", "", .data$chr),
                  X_cum_end = .data$length +.data$ X_cum,
                  text_pos = (.data$length + .data$X_cum*2)/2,
                  Y_end = numofcell-0.8)


  # Fixed color template: Copy number more than 5 as same color
  color_tem <- generate_dynamic_colormap(data_matrix = mat_input$CN)


  # figure legend
  lgd_labels <- sort(unique(mat_input$CN))
  #lgd_labels[which(lgd_labels==5)] <- ">= 5"

  height = dplyr::case_when(
    numofcell > 60 ~ numofcell,
    numofcell <= 60 & numofcell > 10 ~ 50,
    numofcell <= 10  ~ 30
  )

  PlotCN_heatmap <-
    ggplot2::ggplot(mat_input) +
    ggplot2::geom_segment(ggplot2::aes(x = .data$X_cum_start, xend = .data$X_cum_end, y = .data$Y_cum_start, yend = .data$Y_cum_end,
                     color = factor(.data$CN)), linewidth = 2, show.legend = TRUE) +
    ggplot2::scale_color_manual(name="Copy Number", labels=lgd_labels, values=color_tem)+
    ggplot2::geom_segment(data = Chr_tmp, ggplot2::aes(x = .data$X_cum_end , y = -0.2, xend = .data$X_cum_end, yend = .data$Y_end)) +
    ggplot2::geom_text(data = Chr_tmp, ggplot2::aes(x = .data$text_pos, y = -0.05 * height, label = .data$chr), size = 10) +
    ggplot2::labs(x = "Chromosome", y = "Cell ID", color = "Copy Number") +
    ggplot2::xlim(c(0, max(Chr_tmp$X_cum_end)+1)) +
    ggplot2::ylim(c(-0.05 *height, numofcell)) +
    ggplot2::theme(plot.margin = ggplot2::margin(3, 2, 3, 1, "cm"),
                   plot.background = ggplot2::element_rect(fill = "white"),
                   panel.background = ggplot2::element_rect(fill = "white")) +
    ggplot2::theme_void()


  return(PlotCN_heatmap)
}

# 6.3.1: Transform each cell CN-segment information
#' Transform copy number segment information for each cell
#'
#' This function transforms copy number segment information for a single cell,
#' organizing continuous regions with the same copy number into unified segments.
#' It uses chromosome length and cumulative length data to standardize genomic positions.
#'
#' @param data A `GRanges` object representing copy number segments for a single cell.
#' @param index A numeric value indicating the order or ID of the cell. Used for labeling or indexing the output.
#' @param CN_chr_template A data frame containing chromosome information, including chromosome name, chromosome length, and cumulative length.
#'
#' @return A data frame recording copy number segments for the cell, with continuous regions of the same copy number merged.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @keywords internal
#'
segment_transform <- function(data, index, CN_chr_template) {
  as.data.frame(data$bins) %>%
    dplyr::mutate(segment = cumsum(.data$copy.number != stats::lag(.data$copy.number, default = gdata::first(.data$copy.number)))) %>%
    dplyr::group_by(.data$seqnames, .data$segment, .data$copy.number) %>%
    dplyr::summarize(
      start = gdata::first(.data$start),
      end =  gdata::last(.data$end),
      width = (.data$end) - (.data$start) + 1,
      strand = gdata::first(.data$strand),
      .groups = "drop"
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-c(.data$segment, .data$strand)) %>%
    dplyr::rename("chr" = "seqnames", "CN" = "copy.number") %>%
    dplyr::left_join(CN_chr_template, by = "chr") %>%
    dplyr::mutate(
      X_cum_start = .data$start + .data$X_cum,
      X_cum_end = .data$end + .data$X_cum,
      Y_cum_start = index - 1,
      Y_cum_end = index - 1
    )
}

# 6.4: Totalcluster_pdf() function for creating pdf in total clusters by CN matrix
#' Generate copy number profiles for scDNA-seq cell clustering results in PDF file
#'
#' This function outputs copy number profiles from single-cell DNA sequencing (scDNA-seq) clustering results, saving the visualization as a PDF file.
#'
#' @param Input A named list where each element is a `GRanges` object representing a single cell.
#' @param Template A table recorded the clustering, pqArm clustering, re-clustering, and subclone clustering step result for each cell.
#' This table recorded the clustering history in each step.
#' @param FILEname A character string specifying the name of the output PDF file.
#' @param FILEpath A character string specifying the file path where the output PDF will be saved.
#' @param pqArm_file In-build cytoband template for selection: `hg38`, `hg19`, `mm10`, `mm39`.
#' Or a filepath of a table for cytoband information seen on Giemsa-stained chromosomes.
#' It should include the following columns:
#'   - `chrom`: Reference sequence chromosome or scaffold.
#'   - `chromStart`: Start position in genoSeq.
#'   - `chromEnd`: End position in genoSeq.
#'   - `name`: Name of cytogenetic band.
#'   - `gieStain`: Giemsa stain results.
#' @param cellnum_name A character string specifying the column in `Template` that records the number of cells in each cluster.
#' @param cellcutoff A numeric value defining the minimum number of cells required for a cluster to be included.
#' @param cluster_name A character string specifying the column in `Template` that records the clustering result to be visualized.
#' @param sexchromosome_plot A logical value. If `TRUE`, the output plot includes copy number information for sex chromosomes. Defaults to `FALSE`.
#'
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @return A copy number profile visualization saved as a PDF file.
#' @export
#'
#' @examples
#' \dontrun{
#' file_path <- system.file("extdata", "example_data.rds", package = "cnvTree")
#' Example_data <- changeFormat(file = file_path, core = 4)
#' clustering_results <- system.file("extdata", "template_example_data.rds", package = "cnvTree")
#' clustering_results <- readRDS(file = clustering_results)
#'
#' Totalcluster_pdf(
#'   Input = Example_data,
#'   Template = clustering_results,
#'   FILEname = "cnvTree.scDNAseq_Grouping_fig.pdf",
#'   pqArm_file = "hg38",
#'   cellnum_name = "Subclone_cellnum",
#'   cellcutoff = 10,
#'   cluster_name = "Subclone",
#'   sexchromosome_plot = FALSE
#' )
#' }
#'
Totalcluster_pdf <- function(Input, Template, FILEname, FILEpath, pqArm_file, cellnum_name, cellcutoff, cluster_name, sexchromosome_plot=FALSE){
  Cluster_No <- Totalcluster_Cluster_No(Template = Template, cellnum_name = cellnum_name, cellnum = cellcutoff, cluster_name = cluster_name)

  Fig_seq <- list()
  height_ratio <- NULL
  cellnum_list <- NULL
  count = 0
  for (k in Cluster_No){
    cat("Making ", FILEname, ":", "heatmap of", cluster_name, k, "\n")
    SS <- Totalcluster_SS(Template = Template, cluster_name = cluster_name, k = k)
    Cell_num <- length(SS)
    cellnum_list <- c(cellnum_list, Cell_num)


    # height_size <- Cell_num*10
    if(Cell_num<=50){
      height_size <- 300
      height_ratio <- c(height_ratio, 3)
    } else if (Cell_num>50 && Cell_num<=100){
      height_size <- 500
      height_ratio <- c(height_ratio, 5)
    } else {
      height_size <- 1000
      height_ratio <- c(height_ratio, 8)
    }

    count = count + 1
    Fig_seq[[count]] <- GenomeHeatmap(Input = Input, cellID = SS, pqArm_file = pqArm_file, sexchromosome = sexchromosome_plot)

  }
  Label <- paste0(Cluster_No, " (n = ", cellnum_list, ")")


  combine_plot <- ggpubr::ggarrange(plotlist = Fig_seq,
                                    ncol = 1,
                                    labels = Label,
                                    common.legend = FALSE,
                                    legend = "right",
                                    hjust  = 0.8,
                                    align = "v",
                                    font.label = list(size = 45, face = "bold", color ="black"),
                                    heights = c(height_ratio))+
    ggplot2::theme(plot.margin = ggplot2::margin(2,2,2,8, "cm"))



  ggplot2::ggsave(combine_plot,
         filename = paste0(FILEpath, FILEname),
         height = (300*sum(height_ratio) + 800*length(height_ratio)),
         width = 13000,
         units = "px",
         limitsize = FALSE)

}


# 6.4.1: Totalcluster_Cluster_No() function for transfer the number of clusters in the data
#' Identify clusters meeting the thershold by number of cell
#'
#' This function filters clusters based on the number of cells, returning a list of clusters that meet the specified threshold.
#'
#' @param Template  A data frame recording the pqArm clustering, re-clustering, and subclone clustering results for each cell.
#'   This table tracks the clustering history at each step.
#' @param cellnum_name A character string specifying the column in `Template` that records the number of cells in each cluster.
#' @param cellnum A numeric value defining the minimum number of cells required for a cluster to be included.
#' @param cluster_name A character string specifying the column in `Template` that records the clustering result to be filtered.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @return A character vector containing cluster names that meet the specified cell number threshold.
#'
#' @keywords internal
#'
Totalcluster_Cluster_No <- function(Template, cellnum_name, cellnum, cluster_name){
  Cluster_No <- Template %>%
    dplyr::filter(.data[[cellnum_name]] >= cellnum) %>%
    dplyr::pull(.data[[cluster_name]])
  Cluster_No <- sort(unique(Cluster_No))

  return(Cluster_No)
}

# 6.4.2: Totalcluster_SS() function for transfer each cluster of cells in the data
#' Retrieve cell IDs from a specific cluster
#'
#' This function returns a list of cell IDs belonging to a designated cluster based on the clustering results recorded in the `Template`.
#'
#' @param Template A data frame recording the pqArm clustering, re-clustering, and subclone clustering results for each cell.
#'   This table tracks the clustering history at each step.
#' @param cluster_name A character string specifying the column in `Template` that records the clustering result to be queried.
#' @param k An integer specifying the designated cluster for which cell IDs should be retrieved.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @return A character vector containing the cell IDs that belong to the designated cluster.
#'
#' @keywords internal
#'
Totalcluster_SS <- function(Template, cluster_name, k){
  SS <- Template %>%
    dplyr::filter(.data[[cluster_name]] == k) %>%
    dplyr::pull(.data$cellID)

  return(SS)
}


# 6.5: writeOutput() function for outputting cluster results in .txt
#' Export data as a TXT file
#'
#' This function writes a data table to a `.txt` file, saving it to a specified path.
#'
#' @param data A data frame or matrix to be exported as a `.txt` file.
#' @param filename A character string specifying the name of the output `.txt` file.
#' @param path A character string specifying the directory where the file will be saved.
#'
#' @return The function writes a `.txt` file and returns the file path as a character string.
#' @export
#'
writeOutput <- function(data, filename, path){
  FILEpath <- paste0(path, filename, ".txt")
  utils::write.table(data, file = FILEpath, row.names = FALSE, col.names = TRUE)
}


# 6.6: scDNA.superimpose() function for superimpose between definedCNVs and DNA segment info
#' Superimpose defined CNVs onto scDNA-seq copy number results
#'
#' This function overlays high-confidence copy number variations (CNVs) onto
#' single-cell DNA sequencing (scDNA-seq) clustering results to determine
#' whether each cluster contains the corresponding CNVs.
#'
#' @param Template A list containing two data frames:
#'   - `final_cluster_output`: Records the clustering history from the pqArm, re-clustering, and subclone clustering steps.
#'   - `Subclone_CN`: Records each subclone's unique chromosome segment template and its copy number.
#'     It includes the following columns:
#'     - `chr`: Chromosome name (chr1, chr2, ...).
#'     - `start`: Start position of the segment.
#'     - `end`: End position of the segment.
#'     - `region`: The defined region index
#'     - `Subclone`: Subclone identifier.
#'     - `CN`: Copy number value of the segment.
#' @param DefinedCNVs A data frame containing high-confidence CNV regions, with the following columns:
#'   - `chr`: Chromosome name (chr1, chr2, ...).
#'   - `CNV_region`: The index of CNV regions.
#'   - `CN`: Copy number state, categorized as either "amp" (Amplification) or "del" (Deletion).
#'   - `CNV_start`: Start position of the CNV region.
#'   - `CNV_end`: End position of the CNV region.
#'   - `first_band`: Cytoband label of the first affected band.
#'   - `last_band`: Cytoband label of the last affected band.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @return The function returns the updated `Template` list with an additional table:
#'   - `superimpose`: A data frame where `Template$Subclone_CN` overlaps with `DefinedCNVs`,
#'     calculating whether each cluster contains or lacks the corresponding CNV.
#'
#' @keywords internal
#'
scDNA.superimpose <- function(Template, DefinedCNVs){
  Groups <- unique(Template$Subclone_CN$Subclone)
  Template$Subclone_CN$CNV_state <- Template$Subclone_CN$CN
  # CN only seperate in 3 types: del/neu/amp
  Template$Subclone_CN$CN = dplyr::case_when(
    Template$Subclone_CN$CN <  2 ~ "del",
    Template$Subclone_CN$CN == 2 ~ "neu",
    Template$Subclone_CN$CN >  2 ~ "amp")

  superimpose <- NULL
  for(groups in 1:length(Groups)){
    intersection <- NULL
    # check defined CNVs in each group
    CNVs <- Template$Subclone_CN %>%
      dplyr::filter(.data$Subclone %in% c(Groups[groups]),
                    .data$chr %in% c(unique(DefinedCNVs$chr)),
                    .data$CN %in% c(DefinedCNVs$CN))

    intersection <- merge(DefinedCNVs, CNVs, by = c("chr", "CN"))
    intersection <- intersection %>%
      dplyr::mutate(
        F_start = dplyr::case_when(
          .data$start < .data$CNV_start ~ 0,
          .data$start >= .data$CNV_start & .data$start <= .data$CNV_end ~ 1,
          .data$start > .data$CNV_end ~ 2),
        F_end = dplyr::case_when(
          .data$end < .data$CNV_start ~ 0,
          .data$end >= .data$CNV_start & .data$end <= .data$CNV_end ~ 1,
          .data$end > .data$CNV_end ~ 2),
        seg = paste0(.data$F_start, .data$F_end),
        final_start = dplyr::case_when(
          .data$seg %in% c("00", "22") ~ NA,
          .data$seg %in% c("01", "02") ~ .data$CNV_start,
          .data$seg %in% c("11", "12") ~ .data$start),
        final_end = dplyr::case_when(
          .data$seg %in% c("00", "22") ~ NA,
          .data$seg %in% c("01", "11") ~ .data$end,
          .data$seg %in% c("02", "12") ~ .data$CNV_end)) %>%
      dplyr::filter(!.data$seg %in% c("00", "22")) %>%
      dplyr::mutate(cnv_range = .data$final_end - .data$final_start + 1) %>%
      dplyr::group_by(.data$CNV_region) %>%
      dplyr::summarise(cnv_range = sum(.data$cnv_range)) %>%
      as.data.frame()

    superimpose <- dplyr::left_join(DefinedCNVs, intersection, by = "CNV_region") %>%
      dplyr::mutate(Subclone = Groups[groups],
                    CNV_range = .data$CNV_end - .data$CNV_start + 1,
                    cnv_range = ifelse(is.na(.data$cnv_range)==TRUE, 0, .data$cnv_range),
                    cnv_ratio = .data$cnv_range / .data$CNV_range,
                    final_cnv = ifelse(.data$cnv_ratio >= 0.5, 1, 0)) %>%
      rbind(.data, superimpose)
  }

  Template$superimpose <- superimpose


  return(Template)
}


# 6.7: scDNA.clustering() function for receiving DNA clustering output in superimpose range
#' Generate superimposition results between scDNA clustering and defined CNVs
#'
#' This function calculates the overlap between single-cell DNA sequencing (scDNA-seq) clustering results
#' and high-confidence defined copy number variations (CNVs), indicating whether each cluster contains the corresponding CNVs.
#'
#' @param Template A data frame where `Template$Subclone_CN` records each subclone's unique chromosome segment template and its copy number.
#'     It includes the following columns:
#'     - `chr`: Chromosome name (chr1, chr2, ...).
#'     - `start`: Start position of the segment.
#'     - `end`: End position of the segment.
#'     - `region`: The defined region index
#'     - `Subclone`: Subclone identifier.
#'     - `CN`: Copy number value of the segment.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @return The function updates `Template` by adding a new matrix, `Template$DNA_cluster`, which includes:
#'   - `DNA_cluster`: The index for final clustering result.
#'   - `DNA_Cellnum`: The number of cells in each cluster.
#'   - `DefinedCNVs`: A binary matrix (1 or 0) indicating whether each cluster contains the corresponding CNV.
#'
#' @keywords internal
#'
scDNA.clustering <- function(Template){
  cnv_region <- unique(Template$superimpose$CNV_region)
  cnv_region <- paste0("CNV", cnv_region)
  Subclone_ss <- unique(Template$superimpose$Subclone)

  cnv_matrix <- matrix(nrow = length(Subclone_ss), ncol = length(cnv_region)+2)
  for(subclone in 1:length(Subclone_ss)){
    Cell_num <- Template$final_cluster_output %>%
      dplyr::filter(.data$Subclone %in% Subclone_ss[subclone]) %>%
      nrow()
    CNVs <- Template$superimpose %>%
      dplyr::filter(.data$Subclone %in% Subclone_ss[subclone]) %>%
      dplyr::select(.data$final_cnv) %>%
      dplyr::pull()

    cnv_matrix[subclone, ] <-c(CNVs, Subclone_ss[subclone], Cell_num[1])

  }

  colnames(cnv_matrix) <- c(cnv_region, "DNA_cluster", "DNA_Cellnum")
  rownames(cnv_matrix) <- seq_len(nrow(cnv_matrix))

  Template$DNA_cluster <- cnv_matrix


  return(Template)
}

# 6.8: cluster w/ or w/o CNV pattern plot
#' Generate a heatmap of high-confidence CNV patterns with Dendrogram
#'
#' This function creates a heatmap displaying high-confidence copy number variation (CNV) patterns across clusters,
#' with hierarchical clustering represented by a dendrogram.
#'
#' @param Input A named list where each element is a `GRanges` object representing a single cell.
#' @param final_cluster A data frame containing high-confidence CNV regions, with the following columns:
#'   - `chr`: Chromosome name (chr1, chr2, ...).
#'   - `CNV_region`: The index of CNV regions.
#'   - `CN`: Copy number state, categorized as either "amp" (Amplification) or "del" (Deletion).
#'   - `CNV_start`: Start position of the CNV region.
#'   - `CNV_end`: End position of the CNV region.
#'   - `first_band`: Cytoband label of the first affected band.
#'   - `last_band`: Cytoband label of the last affected band.
#' @param cellcutoff A numeric value defining the minimum number of cells required for a cluster to be included.
#' @param pqArm_file In-build cytoband template for selection: `hg38`, `hg19`, `mm10`, `mm39`.
#' Or a filepath of a table for cytoband information seen on Giemsa-stained chromosomes.
#' It should include the following columns:
#'   - `chrom`: Reference sequence chromosome or scaffold.
#'   - `chromStart`: Start position in genoSeq.
#'   - `chromEnd`: End position in genoSeq.
#'   - `name`: Name of cytogenetic band.
#'   - `gieStain`: Giemsa stain results.
#' @param FILEpath A character string specifying the directory where the output file will be saved.
#' @param FILEname A character string specifying the name of the output `.png` file.
#' @param sexchromosome A logical value. If `TRUE`, the heatmap includes copy number information for sex chromosomes.
#' @param smoothing A logical value. If `TRUE`, the heatmap applies smoothing over a 10⁶ bp range in chromosome copy number data.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @return A PNG file containing a heatmap of defined CNVs across clusters.
#' @export
#'
#' @examples
#' \dontrun{
#' file_path <- system.file("extdata", "example_data.rds", package = "cnvTree")
#' Example_data <- changeFormat(file = file_path, core = 4)
#' final_cluster <- system.file("extdata", "template_example_data.rds", package = "cnvTree")
#' final_cluster <- readRDS(file = final_cluster)
#'
#' scDNA_CNVpattern(
#'   Input = Example_data,
#'   final_cluster = final_cluster,
#'   cellcutoff = 5,
#'   pqArm_file = "hg38",
#'   FILEpath = getwd(),
#'   FILEname = "cnvTree.scDNAseq_heatmap.png",
#'   sexchromosome=FALSE,
#'   smoothing = TRUE
#' )
#' }
#'
scDNA_CNVpattern <- function(Input, final_cluster, cellcutoff, pqArm_file, FILEpath, FILEname, sexchromosome=FALSE, smoothing=TRUE){
  if(smoothing==TRUE){
    # bins template
    bins_window = Input[[1]]$bins
    bins_window$copy.number <- NULL

    # chromosome seperated in fixed-bin size
    chr_df = circlize::read.chromInfo()$df
    chr_df = chr_df[chr_df$chr %in% c(paste0("chr", 1:22), "chrX", "chrY"), ] # with sex chromosome
    chr_gr = GenomicRanges::GRanges(seqnames = chr_df[, 1], ranges = IRanges::IRanges(chr_df[, 2]+1, chr_df[, 3]+1))
    chr_window = EnrichedHeatmap::makeWindows(chr_gr, w = 1e6, short.keep = T)

    mtch = as.data.frame(IRanges::findOverlaps(chr_window, bins_window, type = "any"))

    # 先平均出每個bin 的copy number ，得到average sequence
    final_cluster <- final_cluster %>% dplyr::filter(.data$Subclone_cellnum >= cellcutoff)
    subclone_no <- unique(final_cluster$Subclone)
    num_mat <- NULL
    for(i in 1:length(subclone_no)){
      cellID <- final_cluster %>% dplyr::filter(.data$Subclone %in% subclone_no[i]) %>% dplyr::pull(.data$cellID)
      Subclone_CN <- CN_seq(input = Input, Template = cellID)
      averageCN <- round(apply(Subclone_CN, 1, mean), digits = 0)  # mean() for cluster of cells

      smooth_CN = rep(2, length(chr_window))
      recalculateCN = mtch %>%
        dplyr::mutate(copy_number = averageCN[.data$subjectHits]) %>%
        dplyr::group_by(.data$queryHits) %>%
        dplyr::summarise(mean_copy_number = round(mean(.data$copy_number, na.rm = TRUE), digits = 0))

      smooth_CN[as.numeric(recalculateCN$queryHits)] <- recalculateCN$mean_copy_number

      if(is.null(num_mat)){
        num_mat <- smooth_CN
      } else {
        num_mat <- cbind(num_mat, smooth_CN)
      }
    }
    colnames(num_mat) <- subclone_no

    chr_sort <- as.data.frame(chr_window) %>%
      dplyr::select(.data$seqnames, .data$start)
    if(sexchromosome==TRUE) {
      # heatmap annotation labels
      chr <- as.character(sort(GenomicRanges::seqnames(chr_window)))
      chr <- factor(chr, levels = levels(GenomicRanges::seqnames(chr_window)))
    }else if(sexchromosome==FALSE){

      num_mat <- cbind(chr_sort, num_mat) %>%
        dplyr::filter(!.data$seqnames %in% c("chrX", "chrY")) %>%
        as.matrix()
      num_mat <- num_mat[, -c(1, 2)]
      num_mat <- apply(num_mat, c(1, 2), as.numeric)

      # heatmap annotation labels
      chr <- as.character(sort(GenomicRanges::seqnames(chr_window)))
      chr <- chr[!(chr %in% c("chrX", "chrY"))]
      chr <- factor(chr, levels = dplyr::setdiff(levels(GenomicRanges::seqnames(chr_window)), c("chrX", "chrY")))
    }
  } else {
    ### without smoothing step ###
    # 先平均出每個bin 的copy number ，得到average sequence
    final_cluster <- final_cluster %>% dplyr::filter(.data$Subclone_cellnum >= cellcutoff)
    subclone_no <- unique(final_cluster$Subclone)
    num_mat <- NULL
    for(i in 1:length(subclone_no)){
      cellID <- final_cluster %>% dplyr::filter(.data$Subclone %in% subclone_no[i]) %>% dplyr::pull(.data$cellID)
      Subclone_CN <- CN_seq(input = Input, Template = cellID)
      num_mat <- cbind(num_mat, round(apply(Subclone_CN, 1, mean), digits = 0))
    }

    chr_window = Input[[1]]$bins
    chr_window$copy.number <- NULL

    chr_sort <- chr_window %>%
      as.data.frame() %>%
      dplyr::select(.data$seqnames, .data$start) %>%
      dplyr::arrange(.data$seqnames, .data$start)
    if(sexchromosome==TRUE){
      num_mat <- cbind(chr_sort, num_mat) %>%
        dplyr::arrange(.data$seqnames, .data$start) %>%
        as.matrix()

      # heatmap annotation labels
      chr <- as.character(sort(GenomicRanges::seqnames(chr_window)))
      chr <- factor(chr, levels = levels(GenomicRanges::seqnames(chr_window)))
    }else if(sexchromosome==FALSE){
      num_mat <- cbind(chr_sort, num_mat) %>%
        dplyr::arrange(.data$seqnames, .data$start) %>%
        dplyr::filter(!.data$seqnames %in% c("chrX", "chrY")) %>%
        as.matrix()

      # heatmap annotation labels
      chr <- as.character(sort(GenomicRanges::seqnames(chr_window)))
      chr <- chr[!(chr %in% c("chrX", "chrY"))]
      chr <- factor(chr, levels = dplyr::setdiff(levels(GenomicRanges::seqnames(chr_window)), c("chrX", "chrY")))
    }

    num_mat <- num_mat[, -c(1, 2)]
    num_mat <- apply(num_mat, c(1, 2), as.numeric)
  }


  # heatmap annotation labels
  chr_level <- unique(sub("^chr", "", chr))
  dynamic_colors <- generate_dynamic_colormap(data_matrix = num_mat)

  grDevices::png(filename = paste0(FILEpath, FILEname),
      width = 2800,
      height = (length(subclone_no))*150)

  legend_nrow <- min(length(dynamic_colors), 20)
  heatmap_legend <- ComplexHeatmap::Legend(
    labels = names(dynamic_colors),  # 根據您的顏色名稱替換
    legend_gp = grid::gpar(fill = dynamic_colors),
    title = "Copy number",
    nrow = legend_nrow  # 每排顯示10個
  )

  Oncoscan <- ComplexHeatmap::Heatmap(t(num_mat), name = "Copy number", col = dynamic_colors,

                                      column_split = chr,
                                      cluster_columns = FALSE,
                                      cluster_rows = TRUE,
                                      show_row_dend = TRUE,
                                      show_row_names = FALSE,
                                      use_raster = TRUE,
                                      show_heatmap_legend = FALSE,
                                      row_split = subclone_no,

                                      cluster_row_slices = TRUE,

                                      row_dend_width = ggplot2::unit(5, "cm"),
                                      row_dend_gp = grid::gpar(lwd = 2, col = "black"),
                                      row_title = "scDNA clusters",
                                      row_title_gp = grid::gpar(fontsize = 40, fontface = "bold"),
                                      left_annotation = ComplexHeatmap::rowAnnotation(
                                        subgroup = ComplexHeatmap::anno_text(
                                          subclone_no,
                                          rot = 0,
                                          gp = grid::gpar(fontsize = 30))),
                                      column_title = chr_level,
                                      column_title_side = "bottom",
                                      column_title_gp = grid::gpar(fontsize = 25),
                                      border = TRUE,
                                      column_gap = ggplot2::unit(0, "points")
                                      # row_gap = unit(0, "points")
  )

  ComplexHeatmap::draw(Oncoscan,
                       annotation_legend_side = "right",
                       annotation_legend_list = list(heatmap_legend),
                       padding = ggplot2::unit(c(5, 2, 2, 2), "cm"))

  grDevices::dev.off()

}


# 6.8.2: color template for heatmaps (Customize color + RColorBrewer template)
#' Assign colors to copy number states
#'
#' This function maps numerical copy number states to a corresponding color scheme for visualization purposes.
#'
#' @param data_matrix A numeric matrix representing copy number states,
#'   where rows correspond to genomic regions and columns correspond to samples or clusters.
#'
#' @return A matrix of the same dimensions as `data_matrix`, with each numeric copy number
#'   state replaced by a corresponding color code.
#'
#' @keywords internal
#'
generate_dynamic_colormap <- function(data_matrix) {
  unique_values <- sort(unique(as.vector(data_matrix)))

  # 固定的顏色對應表，對應數值 0~5
  predefined_colors <- c("#D0CECE", "#8165A3", "#9BBB59", "#C0504D", "#ED7D31", "#FFC000")
  names(predefined_colors) <- 0:5  # 為 0~5 數值建立顏色對應表

  # 定義漸層調色板，用於處理超過 5 的數值
  gradient_palette <- grDevices::colorRampPalette(c("#c51b7d", "#130303"))

  # 切分數值：0 到 5 和 大於 5
  values_below_6 <- unique_values[unique_values <= 5]
  values_above_5 <- unique_values[unique_values > 5]

  # 分配顏色：0~5 的數值使用對應的顏色
  below_6_colors <- predefined_colors[as.character(values_below_6)]

  # 大於 5 的數值使用漸層顏色
  above_5_colors <- gradient_palette(length(values_above_5))

  # 合併顏色映射
  all_colors <- c(below_6_colors, above_5_colors)
  color_mapping <- stats::setNames(all_colors, c(values_below_6, values_above_5))

  return(color_mapping)
}
