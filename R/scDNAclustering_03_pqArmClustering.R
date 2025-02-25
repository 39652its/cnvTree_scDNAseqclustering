# 3.1: CN_template() build the CN_seq template, stands the range for each row in CN_seq
#' Generate a copy number segment template based on chromosomal arms
#'
#' This function constructs a template for copy number segmentation
#' using cytoband information from Giemsa-stained chromosomes.
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
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @return A data frame containing the defined genomic ranges for p/q arms across all chromosomes.
#' @keywords internal
#'
CN_template <- function(input, pqArm_file){
  CN_tem <- data.frame(GenomicRanges::seqnames(input[[1]]$bins), IRanges::ranges(input[[1]]$bins))
  CN_tem <- CN_tem %>% stats::setNames(c("chr", "start", "end", "width"))

  # Add p/q arm information on the template
  # 非p及q
  pqArm_range <- pqArm_file.remake(FILE = pqArm_file)
  pqArm_range <- pqArm_file.pq(Template = pqArm_range) %>% dplyr::filter(.data$arm == "p")

  CN_tem_pq <- NULL
  for (i in 1:nrow(pqArm_range)){
    CN_tem_pq <- CN_tem %>%
      dplyr::filter(
        .data$chr %in% pqArm_range$chr[i],
        .data$start >= pqArm_range$start[i],
        .data$end <= pqArm_range$end[i]) %>%
      dplyr::mutate(arm = pqArm_range$arm[i]) %>%
      rbind(CN_tem_pq)
  }

  CN_tem_pq <- dplyr::left_join(CN_tem, CN_tem_pq, by = c("chr", "start", "end", "width"))
  CN_tem_pq <- CN_tem_pq %>%
    tidyr::replace_na(list(arm = "q")) %>%
    dplyr::arrange(.data$chr, .data$start)

  return(CN_tem_pq)
}


# 3.1.1: pqArm_file.remake() remake file format from UCSC
#' Extract chromosome arm ranges from UCSC cytoband data
#'
#' This function processes a UCSC cytoband file to extract the ranges of chromosome long
#' and short arms, excluding the centromere regions.
#'
#' @param FILE A character string specifying the file path to the UCSC cytoband file.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @return A table containing the ranges of chromosomal arms, excluding centromeric regions.
#'   It includes cytobands of type `acen` and `gvar`, along with the ranges of the preceding and following cytobands.
#'
#' @keywords internal
#'
pqArm_file.remake <- function(FILE){
  if (FILE=="hg38" | FILE=="hg19" | FILE=="mm10" | FILE=="mm39"){
    filename <- paste0(FILE, "_cytoBand.txt.gz")
    FILE <- system.file("extdata", filename, package = "cnvTree")
  }

  x <- utils::read.table(gzfile(FILE),sep="\t", col.names = c("chr", "start", "end", "name","gieStain"))
  x <- x %>% dplyr::filter(!grepl("_", .data$chr))

  # set chr levels
  vec <- unique(x$chr)
  nums <- as.numeric(gsub("chr", "", vec)[grepl("\\d", vec)])
  nums <- paste0("chr", nums[order(nums)])
  Levels <- c(nums, vec[!grepl("\\d", vec)])


  x <- x %>%
    dplyr::mutate(chr = factor(.data$chr, levels = Levels),
                  start = .data$start + 1,
                  arm  = substring(.data$name, 1, 1),
                  arm_category = paste0(.data$chr, .data$arm))
  Sum_x <- x %>%
    dplyr::group_by(.data$arm_category) %>%
    dplyr::slice_head(n = 1) %>%
    as.data.frame()
  Sum_x <- x %>%
    dplyr::group_by(.data$arm_category) %>%
    dplyr::slice_tail(n = 1) %>%
    as.data.frame() %>%
    rbind(Sum_x) %>%
    dplyr::arrange(chr = factor(.data$chr, levels = Levels), .data$start)

  return(Sum_x)
}


# 3.1.2: pqArm_file.pq() remake pqarm template into new format
#' Convert UCSC cytoband data to arm-level chromosomal ranges
#'
#' This function processes cytoband information from the UCSC database and
#' reformats it into a structured table containing p/q arm regions for each chromosome.
#' The output is designed for downstream copy number variation (CNV) analysis.
#'
#' @param Template A character string specifying the file path to the UCSC cytoband data.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @return A data frame with labeled chromosomal p/q arm regions, formatted for CNV analysis.
#'
#' @keywords internal
#'
pqArm_file.pq <- function(Template){
  x <- Template %>%
    stats::setNames(c("chr", "ChromStart", "ChromEnd", "name", "gieStain", "arm", "arm_category")) %>%
    dplyr::group_by(.data$arm_category) %>%
    dplyr::summarise(start = min(.data$ChromStart),
                     end = max(.data$ChromEnd))
  x <- x %>%
    dplyr::mutate(arm = stringr::str_sub(.data$arm_category, -1),
                  chr = stringr::str_sub(.data$arm_category, end = -2)) %>%
    dplyr::select(c("chr", "start", "end", "arm"))

  return(x)
}

# 3.1.3: pqArm_file.cen() remake centromere template into new format (X)
#' Process UCSC cytoband data for "acen" and "gvar" regions
#'
#' This function extracts and formats cytoband information from the UCSC database,
#' specifically for the "acen" (centromeric) and "gvar" (variable heterochromatic)
#' cytoband types. The output is structured for defining masking ranges
#' to exclude copy number variations (CNVs) from downstream analyses.
#'
#' @param FILE A character string specifying the file path to the UCSC cytoband data file.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @return A data frame containing labeled genomic ranges for "acen" and "gvar" cytoband
#' regions, which can be used for masking CNVs in experimental analyses.
#'
#' @keywords internal
pqArm_file.cen <- function(FILE){
  if (FILE=="hg38" | FILE=="hg19" | FILE=="mm10" | FILE=="mm39"){
    filename <- paste0(FILE, "_cytoBand.txt.gz")
    FILE <- system.file("extdata", filename, package = "cnvTree")
  }

  x <- utils::read.table(gzfile(FILE),sep="\t", col.names = c("chr", "ChromStart", "ChromEnd", "name","gieStain"))
  x <- x %>% dplyr::filter(!grepl("_", .data$chr))

  # set chr levels
  vec <- unique(x$chr)
  nums <- as.numeric(gsub("chr", "", vec)[grepl("\\d", vec)])
  nums <- paste0("chr", nums[order(nums)])
  Levels <- c(nums, vec[!grepl("\\d", vec)])

  x <- x %>%
    dplyr::mutate(cen_category = paste0(.data$chr, .data$gieStain)) %>%
    dplyr::group_by(.data$cen_category) %>%
    dplyr::summarise(Start = min(.data$ChromStart),
                     End  = max(.data$ChromEnd)) %>%
    as.data.frame()

  x <- x %>%
    dplyr::filter(.data$cen_category %in% paste0(rep(Levels,each = 2), c("acen","gvar"))) %>%
    dplyr::mutate(cen = stringr::str_sub(.data$cen_category, -4),
                  chr = stringr::str_sub(.data$cen_category, end = -5)) %>%
    dplyr::group_by(.data$chr) %>%
    dplyr::summarise(MaskStart = min(.data$Start),
                     MaskEnd  = max(.data$End)) %>%
    dplyr::arrange(chr = factor(.data$chr, levels = Levels), .data$MaskStart) %>%
    dplyr::select(c("chr", "MaskStart", "MaskEnd"))

  return(x)
}


# 3.2: pqArm_CN() transform CN matrix to Del/Neu/Amp and based on Arm level to smooth the CN
#' Smooth copy number matrix using arm-level ranges
#'
#' This function processes single-cell copy number data by smoothing
#' copy number variations (CNVs) at the arm level. It utilizes clustering results
#' and cytoband information to assign copy number states (Deletion, Neutral, or Amplification)
#' for each chromosomal arm.
#'
#' @param input A named list where each element is a `GRanges` object representing a single cell.
#' @param Cluster_label An integer specifying the cluster for which the copy number smoothing is performed.
#' @param Clustering_output A data frame recording the clustering results for each cell,
#' including the clustering history at each step.
#' @param pqArm_file In-build cytoband template for selection: `hg38`, `hg19`, `mm10`, `mm39`.
#' Or a filepath of a table for cytoband information seen on Giemsa-stained chromosomes.
#' It should include the following columns:
#'   - `chrom`: Reference sequence chromosome or scaffold.
#'   - `chromStart`: Start position in genoSeq.
#'   - `chromEnd`: End position in genoSeq.
#'   - `name`: Name of cytogenetic band.
#'   - `gieStain`: Giemsa stain results.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @return A matrix with smoothed copy number states (Deletion, Neutral, or Amplification)
#' at the arm level across chromosomes.
#'
#' @keywords internal
#'
pqArm_CN <- function(input, Cluster_label, Clustering_output, pqArm_file){
  selected_files <- subset(Clustering_output, Clustering_output$cluster%in%c(Cluster_label))

  CN_matrix_temp <- CN_template(input = input, pqArm_file = pqArm_file)
  CN_matrix <- CN_seq(input = input, Template = selected_files$cell)
  CN_matrix <- pqArm_DelNeuAmp(matrix = CN_matrix)  # Transform CN to 1, 2, 3

  CN_binsLevel <- CN_matrix_temp %>%
    dplyr::group_by(.data$chr, .data$arm) %>%
    dplyr::summarise(Freq = dplyr::n(), .groups = "drop") %>%
    as.data.frame()
  CN_binsLevel$start <- sapply(1:nrow(CN_binsLevel), function(x){
    sum(CN_binsLevel$Freq[1:x-1])+1
  })
  CN_binsLevel$end <- sapply(1:nrow(CN_binsLevel), function(x){
    sum(CN_binsLevel$Freq[1:x])
  })

  Smooth_pqCN <- NULL
  for(i in 1:nrow(CN_binsLevel)){
    pq_CNmatrix <- CN_matrix[CN_binsLevel$start[i]:CN_binsLevel$end[i], ]
    pq_CN <- sapply(1:ncol(pq_CNmatrix), function(x){
      freq <- table(pq_CNmatrix[ , x]) %>%
        as.data.frame() %>%
        dplyr::arrange(dplyr::desc(.data$Freq)) %>%
        dplyr::pull(.data$Var1) %>%
        as.character() %>%
        as.integer()
      freq[1]
    })

    Smooth_pqCN <- Smooth_pqCN %>%
      rbind(pq_CN)
  }

  Smooth_pqCN <- Smooth_pqCN %>%
    as.data.frame() %>%
    stats::setNames(c(colnames(CN_matrix))) %>%
    `rownames<-`(paste0(rep(levels(CN_matrix_temp$chr), each = 2), rep(c("p", "q"), 11)))

  return(Smooth_pqCN)
}


# 3.2.1: pqArm_DelNeuAmp() transform CN matrix to Del/Neu/Amp three types
#' Categorize copy number variations into three types
#'
#' This function classifies copy number variations (CNVs) into three discrete categories:
#' Deletion (CN < 2), Neutral (CN = 2), and Amplification (CN > 2).
#' The input matrix represents copy number data across genomic regions for multiple cells.
#'
#' @param matrix An integer matrix where columns represent individual cells, and rows correspond
#' to fixed-bin size genomic regions across all chromosomes.
#'
#' @return An integer matrix containing only values 0, 1, and 2, representing:
#'   - `0`: Deletion (CN < 2)
#'   - `1`: Neutral (CN = 2)
#'   - `2`: Amplification (CN > 2)
#'
#' @keywords internal
#'
pqArm_DelNeuAmp <- function(matrix){
  new_matrix <- base::matrix(0, nrow(matrix), ncol(matrix))

  new_matrix[matrix < 2] <- 1
  new_matrix[matrix == 2] <- 2
  new_matrix[matrix > 2] <- 3

  base::colnames(new_matrix) <- base::colnames(matrix)


  return(new_matrix)
}


# 3.3: pqArm_clustering() merge each row as vector to seperate the clsuters
#' Cluster cells based on arm-level copy number patterns
#'
#' This function performs clustering on cells using arm-level copy number variations (CNVs).
#' It groups cells into clusters based on chromosomal arm-level CNV profiles, providing
#' a hierarchical clustering history at each step.
#'
#' @param matrix An integer matrix where columns represent individual cells, and rows correspond
#' to arm-level copy number regions across chromosomes.
#' @param Label An integer specifying the cluster to compute in the arm-level CNV analysis.
#'
#' @return A data frame recording the clustering results for each cell, including
#' the clustering history at each step.
#'
#' @return A table recorded the clustering result for each cell. This table recorded the clustering history in each step.
#' @keywords internal
#'
pqArm_clustering <- function(matrix, Label){
  cluster <- sapply(1:ncol(matrix), function(x){
    paste(matrix[ , x], collapse = "_")
  })

  cluster <- tibble::tibble(
    pqArm_pattern = cluster,
    cellID = colnames(matrix),
    cluster = Label
  )

  return(cluster)
}


# 3.4: pqArm_clustering_summary() return pqArm clustering output
#' Summarize pqArm clustering step results
#'
#' This function provides a summary of the pqArm clustering process,
#' detailing the distribution of arm-level copy number variation (CNV) patterns
#' across different clusters.
#'
#' @param matrix A data frame recording the clustering results for each cell,
#' including the clustering history at each step.
#' @param Label An integer specifying the cluster to compute in the arm-level CNV analysis.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @return A data frame summarizing the pqArm clustering step, containing the following columns:
#'   - `pqArm_pattern`: The identified copy number patterns at the arm level.
#'   - `pqArm_pattern_cellnum`: The number of cells associated with each pqArm pattern.
#'   - `cluster`: The assigned cluster for each pattern.
#'   - `pqArm_cluster`: The final cluster grouping based on pqArm patterns.
#'
#' @keywords internal
#'
pqArm_clustering_summary <- function(matrix, Label){
  cluster_table <- table(matrix$pqArm_pattern) %>%
    as.data.frame() %>%
    dplyr::arrange(dplyr::desc(.data$Freq))
  cluster_table$cluster <- Label
  cluster_table$pqArm_cluster <- seq_len(nrow(cluster_table))
  colnames(cluster_table) <- c("pqArm_pattern", "pqArm_pattern_cellnum", "cluster", "pqArm_cluster")

  return(cluster_table)
}



