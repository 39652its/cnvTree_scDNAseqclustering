# 4.1: pqArm_recluster() calculate the similarity between small clusters and >2 cells clusters
#' Compute similarity between pqArm clusters
#'
#' This function calculates the similarity between pqArm clusters using
#' Euclidean distance. The similarity matrix quantifies the differences
#' in arm-level copy number variation (CNV) patterns across clusters.
#'
#' @param pqArm_cluster A data frame recording the pqArm clustering results for each cell, including the clustering history at each step.
#' @param Cluster An integer specifying the cluster for which similarity calculations are performed.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @return A numeric matrix representing the similarity between pqArm clusters.
#' The similarity is computed based on the Euclidean distance between pqArm cluster patterns.
#'
#' @keywords internal
#'
pqArm_recluster <- function(pqArm_cluster, Cluster){
  # pqArm_cluster <- read.xlsx(xlsxFile = FILEpath) ##這裡需要檢查資料的function
  pqArm_cluster <- pqArm_cluster %>%
    dplyr::filter(.data$cluster%in%Cluster)

  # unique pattern output
  Pattern_more10 <- pqArm_cluster %>%
    dplyr::filter(.data$pqArm_cellnum >= 2) %>%
    dplyr::pull(.data$pqArm_pattern) %>%
    unique()

  # check any Pattern_more10 or Pattern_less10 is NULL
  if(length(Pattern_more10)==0){
    return(0)
  } else{
    Pattern_more10 <- pqArm_cluster.pattern(Pattern_more10)
    # Pattern_less10 <- pqArm_cluster.pattern(pattern = Pattern_less10)
  }


  # New_cluster similarity calculation
  New_cluster <- NULL
  N_cluster <- NULL
  for (i in 1:ncol(Pattern_more10)){
    for (j in 1:ncol(Pattern_more10)){
      N_cluster <- c(N_cluster, euclidean(Pattern_more10[ ,i], Pattern_more10[ ,j]))
    }
    New_cluster <- rbind(New_cluster, N_cluster)
    N_cluster <- NULL
  }

  New_cluster <- New_cluster %>%
    as.data.frame() %>%
    `row.names<-`(colnames(Pattern_more10)) %>%
    `colnames<-`(colnames(Pattern_more10))

  return(New_cluster)

}


# 4.1.1: pqArm_cluster.pattern() convert each cluster pattern from vector to sequence
#' Convert pqArm pattern string to vector format
#'
#' This function transforms a pqArm pattern string into a numeric vector,
#' where each element represents the copy number of a chromosome's p or q arm.
#'
#' @param pattern A character string encoding copy number values for each chromosome's
#' p and q arms, separated by underscores ("_").
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @return A numeric vector containing the parsed copy number values
#' for each chromosome's p and q arm.
#'
#' @keywords internal
#'
pqArm_cluster.pattern <- function(pattern){
  P <- pattern %>%
    base::strsplit(.data , split = "_") %>%
    base::data.frame() %>%
    stats::setNames(c(pattern)) %>%
    dplyr::mutate_if(is.character, as.numeric)

  return(P)
}


# 4.1.2: euclidean() euclidean distance calculation
#' Compute Euclidean distance between two vectors
#'
#' This function calculates the Euclidean distance (L₂ norm) between two numeric vectors.
#' The Euclidean distance is computed as:
#' \deqn{\sqrt{\sum (a_i - b_i)^2}}
#'
#' @param a A numeric vector of the same length as \code{b}.
#' @param b A numeric vector of the same length as \code{a}.
#'
#' @return A numeric value representing the Euclidean distance between the two vectors.
#'
#' @keywords internal
#'
euclidean <- function(a, b){
  sqrt(sum((a - b)^2))
}


# 4.2: pqArm_reclustering_dif() output the different ratio in different pqArm at bins-level between two clusters
#' Compute bin-level difference ratios between pqArm clusters
#'
#' This function calculates the bin-level difference ratios between different
#' pqArm clusters. The ratio represents the degree of difference in copy number
#' variations (CNVs) between each cluster and its most similar cluster,
#' considering variations at the p and q arms.
#'
#' @param input A named list where each element is a `GRanges` object representing a single cell.
#' @param pqArm_recluster_sim A numeric matrix recording the similarity between pqArm clusters.
#' @param pqArm_cluster A data frame recording the pqArm clustering history for each cell.
#' @param Cluster An integer specifying the cluster for which the bin-level
#' difference ratio is computed.
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
#' @return A data frame summarizing bin-level difference ratios across clusters.
#' The output includes:
#'   - `cluster"`: The cluster identifier.
#'   - `pqArm`: The differing p or q arm.
#'   - `bin_difference_ratio`: The ratio of bin-level differences between each cluster and its most similar cluster.
#'
#' @keywords internal
#'
pqArm_reclustering_dif <- function(input, pqArm_recluster_sim, pqArm_cluster, Cluster, pqArm_file){
  pqArm_sim <- apply(pqArm_recluster_sim, 2, function(x) min(x[x!=0]) )   # 2: column is the more10 cluster pattern


  pqArm_sim <- as.data.frame(pqArm_sim)

  CN_bins_template <- CN_template(input = input, pqArm_file = pqArm_file)
  chr_pq <- paste0(rep(levels(CN_bins_template$chr), each = 2), rep(c("p", "q"), 11))

  # 將細胞數<10的cluster 與細胞數>2的cluster 相比，找到相似度最高的>2 Cluster 考慮合併
  new_pqArm_cluster <- NULL
  for(i in 1:nrow(pqArm_sim)){
    p <- rownames(pqArm_sim)[i]
    more10_Row <- which(pqArm_recluster_sim[,p] == pqArm_sim[i, 1])
    less10 <- rep(p, times = length(more10_Row))
    new_pqArm_cluster <- cbind(less10, rownames(pqArm_recluster_sim)[more10_Row]) %>%
      rbind(new_pqArm_cluster)
  }
  new_pqArm_cluster <- new_pqArm_cluster %>%
    as.data.frame() %>%
    stats::setNames(c("less10", "more10"))

  # 得到 information about which pqArm is different
  new_pqArm_PQreturn <- NULL
  for(i in 1:nrow(new_pqArm_cluster)){
    p <- c(new_pqArm_cluster$less10[i], new_pqArm_cluster$more10[i])
    PQreturn <- pqArm_return.PQ(pattern = p, PQarm = chr_pq)  # names(new_pqArm_PQreturn): Subclone name, Inside: CellNum>10 cluster pattern
    Times <- length(PQreturn)
    lessmore <- cbind(less10 = rep(new_pqArm_cluster$less10[i], times = Times), more10 = rep(new_pqArm_cluster$more10[i], times = Times))
    new_pqArm_PQreturn <- cbind(lessmore, PQreturn) %>%
      rbind(new_pqArm_PQreturn)%>%
      as.data.frame()
  }


  # which pqArm is different than change into bins-level than check how many bins are different
  # pqArm_cluster <- read.xlsx(xlsxFile = FILEpath)
  pqArm_cluster <- pqArm_cluster %>%
    dplyr::filter(.data$cluster == Cluster)
  CN_matrix <- CN_seq(input = input, Template = pqArm_cluster$cellID)
  CN_matrix$Chr_arm <- paste0(CN_bins_template$chr, CN_bins_template$arm)

  dif_num <- NULL
  dif_ratio <- NULL
  for(i in 1:nrow(new_pqArm_PQreturn)){
    Arm = new_pqArm_PQreturn$PQreturn[i]
    more10_CN <- pqArm_return.Bins(Pattern = new_pqArm_PQreturn$more10[i], which_Arm = Arm, Tem = pqArm_cluster, CN_matrix = CN_matrix)
    less10_CN <- pqArm_return.Bins(Pattern = new_pqArm_PQreturn$less10[i], which_Arm = Arm, Tem = pqArm_cluster, CN_matrix = CN_matrix)

    dif_num <- c(dif_num, length(which(more10_CN != less10_CN))) # Number of bins are different
    dif_ratio <- c(dif_ratio, length(which(more10_CN != less10_CN))/length(more10_CN)) # Ratio in chr are different
  }
  new_pqArm_PQreturn$dif_num <- dif_num
  new_pqArm_PQreturn$dif_ratio <- dif_ratio

  return(new_pqArm_PQreturn)
}


# 4.2.1: pqArm_return.PQ() output the different pqArm between two clusters
#' Identify differences in p/q arm positions between chromosome clusters
#'
#' This function detects differences in the p or q arm positions of chromosomes
#' between two clusters based on copy number variation (CNV) patterns.
#'
#' @param pattern A data frame with two columns representing paired CNV patterns.
#'   Each row corresponds to a specific chromosome region comparison.
#' @param PQarm A character vector specifying the chromosome order list.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @return A character vector indicating the chromosome arms (`p` or `q`)
#'   that show differences between the two clusters.
#'
#' @keywords internal
#'
pqArm_return.PQ <- function(pattern, PQarm){
  pqArm_list <- PQarm

  Pattern_unlist <- pqArm_cluster.pattern(pattern = pattern) %>%
    stats::setNames(c("less10", "more10"))

  pqArm_select <- pqArm_list[which(Pattern_unlist$less10 != Pattern_unlist$more10)]


  return(pqArm_select)
}


# 4.2.2: pqArm_return.Bins() select different pqArm to output the region at bin-level
#' Extract bin-level copy number sequence in differentiated chromosome arms
#'
#' This function retrieves the bin-level copy number sequence for chromosome p/q arms
#' that show differences between two clusters.
#'
#' @param Pattern A character string encoding copy number values for each chromosome p/q arm,
#'   with values concatenated using `_`.
#' @param which_Arm A character vector specifying the chromosome arms (`p` or `q`)
#'   that differ between the two clusters.
#' @param Tem A data frame recording the clustering history for each cell, including
#'   overall clustering results and p/q arm-specific clustering assignments.
#' @param CN_matrix An integer matrix where:
#'   - Columns represent the `cellID`s of the specific cells.
#'   - Rows represent genomic regions, divided into fixed bins.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @return A numeric vector where each element represents the mode copy number value
#'   for a given bin within a cluster of cells.
#'
#' @keywords internal
#'
pqArm_return.Bins <- function(Pattern, which_Arm, Tem, CN_matrix){
  ID <- Tem %>%
    dplyr::filter(.data$pqArm_pattern %in% c(Pattern)) %>%
    dplyr::pull(.data$cellID)

  if (length(ID) == 0) {
    stop("No matching cellID found in copy number matrix.")
  }

  CNmatrix <- CN_matrix %>%
    dplyr::filter(.data$Chr_arm %in% which_Arm) %>%
    dplyr::select(dplyr::all_of(ID))  # 確保 ID 為存在的列名

  CNmatrix <- pqArm_DelNeuAmp(matrix = CNmatrix)  # 只看 Del/Neu/Amp

  CN_bins <- apply(CNmatrix, 1, function(x) {
    freq <- table(as.integer(x))
    sorted_freq <- sort(freq, decreasing = TRUE)
    first_element <- as.integer(names(sorted_freq)[1])
    return(first_element)
  })


  return(CN_bins)
}


# 4.3: pqArm_reclusterBy_ratio_target() filter ratio and merge the clusters if criteria meets
#' Validate pqArm clusters based on difference ratio criteria
#'
#' This function evaluates whether pqArm clusters meet a predefined difference ratio criterion across different chromosomes.
#'
#' @param pqArm_cluster A data frame recording the pqArm-specific clustering results for each cell.
#'   This table tracks the clustering history at each step.
#' @param Cluster An integer specifying the cluster to analyze.
#' @param pqReclsut_sim A data frame containing the similarity scores between the entire cluster and each different ratio at the bin level.
#'   It returns the bin-level difference ratio between each cluster and its most similar cluster for each differing p or q arm.
#' @param difratio_chr A numeric value defining the threshold for acceptable difference ratios across different chromosomes.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @return A list containing character vectors representing clusters that should be grouped together.
#'
#' @keywords internal
#'
pqArm_reclusterBy_ratio_target <- function(pqArm_cluster, Cluster, pqReclsut_sim, difratio_chr){
  # pqArm_cluster <- read.xlsx(xlsxFile = FILEpath)
  pqArm_cluster <- pqArm_cluster %>%
    dplyr::filter(.data$cluster %in% Cluster)

  Chioce <- pqReclsut_sim %>%
    dplyr::group_by(.data$less10) %>%
    dplyr::mutate(less10_times = dplyr::n(),
                  merge_pattern = paste0(.data$less10, "_", .data$more10))
  cellnum <- pqArm_cluster %>%
    dplyr::select(.data$pqArm_pattern, .data$pqArm_cellnum) %>%
    dplyr::filter(.data$pqArm_pattern %in% c(Chioce$more10)) %>%
    dplyr::distinct(.data$pqArm_pattern, .keep_all = TRUE) %>%
    stats::setNames(c("more10", "more10_cellnum"))
  Chioce <- merge(Chioce, cellnum, by = "more10")


  # 多個region 不同的要都符合才能留下
  for (pattern in unique(Chioce$merge_pattern)){
    Selected <- Chioce %>%
      dplyr::filter(.data$merge_pattern %in% c(pattern),
                    .data$dif_ratio > difratio_chr)

    if(nrow(Selected)>0){
      Chioce <- Chioce %>%
        dplyr::filter(!.data$merge_pattern %in% c(pattern))
    } else {
      Chioce <- Chioce
    }
  }

  # 一種less10 最終只能配對到一個more10
  Chioce_Result <- NULL
  #pattern = unique(Chioce$less10)[3]
  for (pattern in unique(Chioce$less10)){
    Selected <- Chioce %>%
      dplyr::select(.data$less10, .data$more10, .data$PQreturn, .data$dif_num, .data$dif_ratio, .data$more10_cellnum) %>%
      dplyr::filter(.data$less10 %in% c(pattern))
    if(length(unique(Selected$more10))>1){
      Selected <- Selected %>%
        dplyr::filter(.data$dif_ratio == min(Selected$dif_ratio))
      Chioce_Result <- Chioce_Result %>%
        rbind(Selected)
    } else{
      Chioce_Result <- Chioce_Result %>%
        rbind(Selected)
    }
  }

  if (is.null(Chioce_Result) == TRUE){
    return(NULL)
  } else {
    new_Chioce <- list()
    count = 0
    while(nrow(Chioce_Result)>0){
      count = count + 1
      pattern <- c(Chioce_Result$less10[1], Chioce_Result$more10[1]) # select start merge cluster
      ss <- Chioce_Result %>%
        dplyr::filter(.data$less10 %in% pattern | .data$more10 %in% pattern)
      pattern_group <- c(unique(ss$less10, ss$more10))
      new_Chioce[[count]] <- pattern_group

      Chioce_Result <- Chioce_Result %>%
        dplyr::filter(!.data$less10 %in% pattern_group & !.data$more10 %in% pattern_group)
    }
  }

  return(new_Chioce)
}


# 4.4: pqArm_recluster_result() reset the content in pqArmCluster_CellID.xlsx and merge final result
#' Generate final results of the Re-Clustering step
#'
#' This function updates the clustering results by incorporating the re-clustering step,
#' refining the pqArm-based clustering assignments.
#'
#' @param pqArm_cluster A data frame recording the pqArm-specific clustering results for each cell.
#'   This table tracks the clustering history at each step.
#' @param pqReclsut_target A list containing character vectors representing clusters that should be grouped together.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @return A data frame containing updated clustering results, including:
#'   - `cluster`: The original clustering assignments.
#'   - `pq_cluster`: The pqArm-based clustering assignments.
#'   - `recluster`: The updated clustering result after the re-clustering step.
#'
#' @keywords internal
#'
pqArm_recluster_result <- function(pqArm_cluster, pqReclsut_target){
  pqArm_cluster$Recluster_pattern <- pqArm_cluster$pqArm_pattern

  for(i in 1:length(pqReclsut_target)){
    pqArm_cluster$Recluster_pattern <- ifelse(pqArm_cluster$Recluster_pattern%in%pqReclsut_target[[i]], i, pqArm_cluster$Recluster_pattern)
  }

  Recluster_summary <- pqArm_reclustering_summary(Data = pqArm_cluster$Recluster_pattern)
  pqArm_cluster <- dplyr::left_join(pqArm_cluster, Recluster_summary , by = "Recluster_pattern") %>%
    dplyr::arrange(dplyr::desc(.data$Recluster_cellnum)) %>%
    dplyr::select(.data$cellID, .data$cluster, .data$pqArm_pattern, .data$pqArm_cellnum, .data$pqArm_cluster,
                  .data$Recluster_pattern, .data$Recluster_cellnum, .data$Recluster_cluster)

  return(pqArm_cluster)
}


# 4.4.1: pqArm_reclustering_summary() create pqArm clustering final results
#' Summarize the results of the Re-Clustering step
#'
#' This function generates a summary table of the re-clustering process,
#' providing an overview of the reclustered patterns, cell counts, and assigned clusters.
#'
#' @param Data A data frame recording the clustering results and pqArm-specific clustering results for each cell.
#'   This table tracks the clustering history at each step.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @return A data frame summarizing the re-clustering step with the following columns:
#'   - `Recluster_pattern`: The identified re-clustering patterns.
#'   - `Recluster_cellnum`: The number of cells assigned to each pattern.
#'   - `Recluster_cluster`: The final cluster assignment after re-clustering.
#'
#' @keywords internal
#'
pqArm_reclustering_summary <- function(Data){
  cluster_table <- table(Data) %>%
    as.data.frame() %>%
    dplyr::arrange(dplyr::desc(.data$Freq))
  cluster_table$Recluster_cluster <- seq_len(nrow(cluster_table))
  colnames(cluster_table) <- c("Recluster_pattern", "Recluster_cellnum", "Recluster_cluster")

  return(cluster_table)
}
