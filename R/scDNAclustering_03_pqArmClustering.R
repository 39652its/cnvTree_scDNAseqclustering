# 3.1: CN_template() build the CN_seq template, stands the range for each row in CN_seq
#' Build the template for copy number segment template
#'
#' @param input A list of cells in GRanges-format object.
#' @param pqArm_file A table for cytoband information seen on Giemsa-stained chromosomes.
#'
#' @return A table recorded ranges of p/q arm on each chromosome.
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#'
CN_template <- function(input, pqArm_file){
  CN_tem <- data.frame(GenomicRanges::seqnames(input[[1]]$bins), IRanges::ranges(input[[1]]$bins))
  CN_tem <- CN_tem %>% dplyr::setNames(c("chr", "start", "end", "width"))

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
#' Remake the cytoband information from UCSC database for "acen" and "gvar" cytoband types to specific format
#'
#' @param FILE A directory direct to UCSC cytoband file.
#'
#' @return A table labeled ranges of "acen" and "gvar" cytoband types for further masking ranges to exclude the CNVs from experiments.
#' @export
#'
#' @examples
pqArm_file.remake <- function(FILE){
  x <- utils::read.table(gzfile(FILE),sep="\t", col.names = c("chr", "start", "end", "name","gieStain"))
  x <- x %>% dplyr::filter(!grepl("_", chr))

  # set chr levels
  vec <- unique(x$chr)
  nums <- as.numeric(gsub("chr", "", vec)[grepl("\\d", vec)])
  nums <- paste0("chr", nums[order(nums)])
  Levels <- c(nums, vec[!grepl("\\d", vec)])


  x <- x %>%
    dplyr::mutate(chr = factor(chr, levels = Levels),
                  start = start + 1,
                  arm  = substring(name, 1, 1),
                  arm_category = paste0(chr, arm))
  Sum_x <- x %>%
    dplyr::group_by(arm_category) %>%
    dplyr::slice_head(n = 1) %>%
    as.data.frame()
  Sum_x <- x %>%
    dplyr::group_by(arm_category) %>%
    dplyr::slice_tail(n = 1) %>%
    as.data.frame() %>%
    rbind(Sum_x) %>%
    dplyr::arrange(chr = factor(chr, levels = Levels), start)

  return(Sum_x)
}


# 3.1.2: pqArm_file.pq() remake pqarm template into new format
#' Remake the cytoband information from UCSC database for arm-level ranges to specific format
#'
#' @param Template A directory direct to UCSC cytoband file.
#'
#' @return A table labeled ranges of p/q arm regions in each chromosome for further CNVs calculation.
#' @export
#'
#' @examples
pqArm_file.pq <- function(Template){
  x <- Template %>%
    dplyr::setNames(c("chr", "ChromStart", "ChromEnd", "name", "gieStain", "arm", "arm_category")) %>%
    dplyr::group_by(arm_category) %>%
    dplyr::summarise(start = min(ChromStart),
                     end = max(ChromEnd))
  x <- x %>%
    dplyr::mutate(arm = stringr::str_sub(arm_category, -1),
                  chr = stringr::str_sub(arm_category, end = -2)) %>%
    dplyr::select(c("chr", "start", "end", "arm"))

  return(x)
}


# 3.2: pqArm_CN() transform CN matrix to Del/Neu/Amp and based on Arm level to smooth the CN
#' Smooth copy number matrix based on arm-level ranges
#'
#' @param input A list of cells in GRanges-format object.
#' @param Cluster_label A integer for the specific cluster to calculate.
#' @param Clustering_output A table recorded the clustering result for each cell. This table recorded the clustering history in each step.
#' @param pqArm_file A table for cytoband information seen on Giemsa-stained chromosomes.
#'
#' @return A matrix in three copy number types (Deletion/Neutral/Amplification) based on arm-level ranges.
#' @export
#'
#' @examples
pqArm_CN <- function(input, Cluster_label, Clustering_output, pqArm_file){
  selected_files <- subset(Clustering_output, Clustering_output$cluster%in%c(Cluster_label))

  CN_matrix_temp <- CN_template(input = input, pqArm_file = pqArm_file)
  CN_matrix <- CN_seq(input = input, Template = selected_files$cell)
  CN_matrix <- pqArm_DelNeuAmp(matrix = CN_matrix)  # Transform CN to 1, 2, 3

  CN_binsLevel <- CN_matrix_temp %>%
    dplyr::group_by(chr, arm) %>%
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
        dplyr::arrange(dplyr::desc(Freq)) %>%
        dplyr::pull(Var1) %>%
        as.character() %>%
        as.integer()
      freq[1]
    })

    Smooth_pqCN <- Smooth_pqCN %>%
      rbind(pq_CN)
  }

  Smooth_pqCN <- Smooth_pqCN %>%
    as.data.frame() %>%
    dplyr::setNames(c(colnames(CN_matrix))) %>%
    `rownames<-`(paste0(rep(levels(CN_matrix_temp$chr), each = 2), rep(c("p", "q"), 11)))

  return(Smooth_pqCN)
}


# 3.2.1: pqArm_DelNeuAmp() transform CN matrix to Del/Neu/Amp three types
#' Based on copy number to separate in three types
#'
#' @param matrix An integer matrix which columns are cells and rows are fixed-bin size region in whole chromosomes.
#'
#' @return A integer matrix only with 0, 1, 2 stands for Deletion(CN<2), Neutral(CN=2), and Amplification(CN>2).
#' @export
#'
#' @examples
pqArm_DelNeuAmp <- function(matrix){
  new_matrix <- base::matrix(0, nrow(matrix), ncol(matrix))

  new_matrix[matrix < 2] <- 1
  new_matrix[matrix == 2] <- 2
  new_matrix[matrix > 2] <- 3

  base::colnames(new_matrix) <- base::colnames(matrix)


  return(new_matrix)
}


# 3.3: pqArm_clustering() merge each row as vector to seperate the clsuters
#' Based on arm-level copy number pattern to cluster the cells
#'
#' @param matrix An integer matrix which columns are cells and rows are arm-level in each chromosome.
#' @param Label An integer for the specific cluster to calculate in arm-level result.
#'
#' @return A table recorded the clustering result for each cell. This table recorded the clustering history in each step.
#' @export
#'
#' @examples
pqArm_clustering <- function(matrix, Label){
  cluster <- sapply(1:ncol(matrix), function(x){
    paste(matrix[ , x], collapse = "_")
  })

  cluster <- cluster %>%
    as.data.frame() %>%
    tibble::tibble(pqArm_pattern = ., cellID = colnames(matrix)) %>%
    dplyr::mutate(cluster = Label)


  return(cluster)
}


# 3.4: pqArm_clustering_summary() return pqArm clustering output
#' Return summary of pqArm clustering step result
#'
#' @param matrix A table recorded the clustering result for each cell. This table recorded the clustering history in each step.
#' @param Label An integer for the specific cluster to calculate in arm-level result.
#'
#' @return A table with columns "pqArm_pattern", "pqArm_pattern_cellnum", "cluster", and "pqArm_cluster" as summary of pqArm clustering step.
#' @export
#'
#' @examples
pqArm_clustering_summary <- function(matrix, Label){
  cluster_table <- table(matrix$pqArm_pattern) %>%
    as.data.frame() %>%
    dplyr::arrange(dplyr::desc(Freq)) %>%
    dplyr::mutate(cluster = Label,
                  pqArm_cluster = c(1:nrow(.))) %>%
    # filter(Freq > 1) %>%
    dplyr::setNames(c("pqArm_pattern", "pqArm_pattern_cellnum", "cluster", "pqArm_cluster"))

  return(cluster_table)
}



