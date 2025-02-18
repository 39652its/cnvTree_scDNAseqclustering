##### 02: Clustering

# 2.0: clusterbyHMM() based on AneuFinder::clusterHMMs(), calcuate distance than hierrachical clustering
#' Hierarchical clustering of cells based on copy number variation
#'
#' This function computes the pairwise distance between cells based on their copy number variations
#' and performs hierarchical clustering.
#'
#' @param input A named list where each element is a `GRanges` object representing a single cell.
#' @param selected A character vector specifying the `cellID`s of the cells to be included in clustering.
#' @param exclude.regions A `GRanges` object specifying genomic regions to exclude from clustering computation.
#'   This is useful for filtering out regions with artifacts.
#'
#' @return A list containing:
#'   - `ordered_indices`: The ordered indices of cells based on hierarchical clustering.

#' @export
#'
#' @examples
#'
#' file_path <- system.file("extdata", "example_data.rds", package = "cnvTree")
#' Example_data <- changeFormat(file = file_path, core = 4)
#' Clustering_Result <- clusterbyHMM(input = Example_data, selected = names(Example_data)[1:10])
#'
#'
clusterbyHMM <- function(input, selected, exclude.regions = NULL){
  hmms <- input[selected]

  ptm <- startTimed("Checking column 'copy.number'  ...")
  hmms2use <- numeric()
  for (i1 in 1:length(hmms)) {
    hmm <- hmms[[i1]]
    if (!is.null(hmm$bins$copy.number)) {
      if (is.null(hmm$ID)) {
        stop("Need ID to continue.")
      }
      hmms2use[hmm$ID] <- i1
    }
  }
  hmms <- hmms[hmms2use]
  endTimed(ptm)

  hc <- NULL

  ptm <- startTimed("Making consensus template ...")
  if (!is.null(hmms[[1]]$bins$copy.number)) {
    constates <- sapply(hmms, function(hmm) {
      hmm$bins$copy.number
    })
  }
  constates[is.na(constates)] <- 0
  vars <- apply(constates, 1, stats::var, na.rm = TRUE)
  endTimed(ptm)

  ptm <- startTimed("Clustering ...")
  if (!is.null(exclude.regions)) {
    ind <- GenomicRanges::findOverlaps(hmms[[1]]$bins, exclude.regions)@from
    constates <- constates[-ind, ]
  }

  # ptm <- startTimed("Distance calculating...")
  Dist <- Rfast::Dist(t(constates), method = "euclidean")
  # endTimed(ptm)

  # dist <- parallelDist::parDist(t(constates),
  #                               method = "euclidean",
  #                               threads = 5) # threads

  Dist_as_dist <- stats::as.dist(Dist)
  # ptm <- startTimed("hierarchical clustering...")
  hc <- stats::hclust(Dist_as_dist)
  endTimed(ptm)

  # message("Reordering ...")
  hmms2use <- hmms2use[hc$order]


  return(list(IDorder = hmms2use, hclust = hc))
}

# 2.1: CutTree_final() get the phylogenetic tree template
#' Construct a phylogenetic tree from copy number variation data
#'
#' This function performs hierarchical clustering on selected cells based on their copy number variations
#' and derives a phylogenetic tree structure.
#'
#' @param input A named list where each element is a `GRanges` object representing a single cell.
#' @param selected A character vector specifying the `cellID`s of the cells to be included in the phylogenetic analysis.
#'
#' @return A data frame with two columns:
#'   - `cellID`: The unique identifier of each cell.
#'   - `cluster`: The assigned cluster label (k = 2).
#'
#' @keywords internal
#'
CutTree_final <- function(input, selected){
  # 分群的原始檔，後面要用他作為基底
  message("Clustering and Data processing ...")

  clust <- clusterbyHMM(input = input, selected = selected)
  Clust_cuttree <- data.frame(cluster = stats::cutree(clust[["hclust"]], k = 2),
                              cell = names(clust$IDorder))

  Clust_cuttree <- Clust_cuttree %>%
    dplyr::mutate(cluster = as.numeric(.data$cluster)) %>%
    dplyr::arrange(dplyr::desc(.data$cluster))


  return(Clust_cuttree)
}

# 2.2: CutTree() Cut the tree repeatably
#' Divide cells into two groups based on copy number variation in a specific cluster
#'
#' This function separates cells into two groups based on copy number variation in a specified cluster label
#' from a provided clustering result table.
#'
#' @param input A named list where each element is a `GRanges` object representing a single cell.
#' @param Template A data frame containing two columns:
#'   - `cellID`: Unique identifier for each cell.
#'   - `cluster`: Cluster assignment for each cell.
#' @param Cluster_label An integer specifying the cluster label used to divide the cells into two groups.
#'
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#'
#' @return A data frame with two columns:
#'   - `cellID`: The unique identifier of each cell.
#'   - `cluster`: The updated cluster assignment (k = 2).
#'
#' @keywords internal
#'
CutTree <- function(input, Template, Cluster_label){
  # selected.files建立
  # Cluster_label: Clust_cuttree$cluster中的分群數字
  selected.files <- NULL
  selected.files <- subset(Template, .data$cluster%in%c(Cluster_label))$cell

  message("Divided the cluster in k=2 ")
  clust <- clusterbyHMM(input = input, selected = selected.files, exclude.regions = NULL)
  Clust_2 <- data.frame(NewCluster = stats::cutree(clust[["hclust"]], k = 2),
                        cell = names(clust$IDorder))


  Template <- merge(Template, Clust_2, by = "cell", all = TRUE)
  Template$NewCluster <- ifelse(is.na(Template$NewCluster)==T, 0, Template$NewCluster)

  Clust_2 <- which(is.na(Template$NewCluster) == F)

  count = max(Template$cluster, na.rm = TRUE)

  Template$cluster <- dplyr::case_when(
    Template$NewCluster == 1 ~ count + 1,
    Template$NewCluster == 2 ~ count + 2,
    TRUE ~ Template$cluster
  )

  Template <- Template[, !colnames(Template) %in% "NewCluster"]

  return(Template)
}

# 2.3: Cluster_num() calculate the number of cells in cluster
#' Count the number of cells in a specific cluster
#'
#' This function calculates the total number of cells that belong to a specified cluster.
#'
#' @param Template A data frame containing two columns:
#'   - `cellID`: Unique identifier for each cell.
#'   - `cluster`: Cluster assignment for each cell.
#' @param Cluster_label An integer specifying the cluster for which the number of cells will be counted.
#'
#' @return An integer representing the number of cells in the specified cluster.
#'
#' @keywords internal
Cluster_num <- function(Template, Cluster_label){
  # cat("Calculating numbers of cell in Cluster", Cluster_label, "...\n")

  num <- data.frame(table(Template$cluster))
  num <- num[which(num$Var1 == Cluster_label), 2]

  return(num)
}

# 2.4: Cluster_sim() calculate the similarity of cells in cluster
#' Compute the similarity of cells within a cluster
#'
#' This function calculates the overall similarity of cells within a specified cluster
#' using a precomputed cell-to-cell similarity matrix.
#'
#' @param Template A data frame containing two columns:
#'   - `cellID`: Unique identifier for each cell.
#'   - `cluster`: Cluster assignment for each cell.
#' @param SimCells A square matrix where each element `[i, j]` represents the similarity
#'   score between `cellID[i]` and `cellID[j]`.
#' @param Cluster_label A character vector containing the `cellID`s of cells that belong to the specified cluster.
#'
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#'
#' @return A numeric value representing the overall similarity of the cells within the specified cluster.
#'
#' @keywords internal
#'
Cluster_sim <- function(Template, SimCells, Cluster_label){
  # cat("Calculating cell similarity in Cluster ",  Cluster_label, " ...\n")

  selected <- Template %>% dplyr::filter(.data$cluster %in% Cluster_label) %>% dplyr::pull(.data$cell)
  selected <- which(SimCells$IDorder %in% selected)

  Similarity <- SimCells$similarity[selected, selected]
  Similarity[is.na(Similarity)] <- 0
  Similarity <- mean(Similarity)


  return(Similarity)
}


#' Compute pairwise similarity between cells
#'
#' This function calculates the similarity values between cells based on their copy number
#' variations at the bin level.
#'
#' @param binsMatrix A numeric matrix where each row represents a genomic bin
#'   and each column represents a cell. The values indicate copy number variations.
#'
#' @return A square numeric matrix where each element `[i, j]` represents the similarity
#'   score between `cell[i]` and `cell[j]`.
#'
#' @keywords internal
#'
Cluster_SimTem <- function(binsMatrix) {
  ptm <- startTimed("Making similarity template ... ")

  totalcells <- ncol(binsMatrix)
  num_bins <- nrow(binsMatrix)
  result <- sapply(1:totalcells, function(i) {
    sapply(i:totalcells, function(j) {
      sum(binsMatrix[, i] == binsMatrix[, j]) / num_bins
    })
  })

  similarity <- matrix(0, nrow = totalcells, ncol = totalcells)
  for (i in 1:totalcells) {
    similarity[i, i:totalcells] <- result[[i]]
    similarity[i:totalcells, i] <- result[[i]]
  }

  endTimed(ptm)

  SIM <- list(IDorder = colnames(binsMatrix),
              similarity = as.matrix(similarity))

  return(SIM)
}
