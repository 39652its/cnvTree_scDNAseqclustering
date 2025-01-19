##### 02: Clustering

# 2.0: clusterbyHMM() based on AneuFinder::clusterHMMs(), calcuate distance than hierrachical clustering
#' Calculate each cell distance by each cell copy number than heirrachical clustering
#'
#' @param input A list of cells in GRanges-format object.
#' @param selected A list of cell IDs for cell clustering.
#' @param exclude.regions A GRanges-class with regions that will be excluded from the computation of the clustering. This can be useful to exclude regions with artifacts.
#'
#' @return A list() with ordered ID indices and the hierarchical clustering.
#' @export
#'
#' @examples
#' file_path <- system.file("extdata", "example_data.rds", package = "cnvTree")
#' Example_data <- changeFormat(file = file_path, core = 4)
#' Clustering_Result <- clusterbyHMM(input = Example_data, selected = names(Example_data)[1:10])
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
  vars <- apply(constates, 1, var, na.rm = TRUE)
  endTimed(ptm)

  ptm <- startTimed("Clustering ...")
  if (!is.null(exclude.regions)) {
    ind <- findOverlaps(hmms[[1]]$bins, exclude.regions)@from
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
#' Calculate phylogenetic tree results
#'
#' @param input A list of cells in GRanges-format object.
#' @param selected A list of cell IDs for cell clustering.
#'
#' @return A table recorded cell IDs and number of cluster(k=2).
#' @export
#'
#' @examples
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
#' Divided in 2 groups based on the specific cluster
#'
#' @param input A list of cells in GRanges-format object.
#' @param Template A table result contains cellIDs and cluster column.
#' @param Cluster_label An integer. The specific cluster to divide in 2 cluster.
#'
#' @return A table recorded cell IDs and number of cluster(k=2).
#' @export
#'
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#'
#' @examples
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
#' Calculate the number of cell in cluster
#'
#' @param Template A table result contains cellIDs and cluster column.
#' @param Cluster_label An integer. The specific cluster need to calculate the number of cells.
#'
#' @return An integer represent the number of cell in specific cluster.
#' @export
#'
#' @examples
#'
Cluster_num <- function(Template, Cluster_label){
  # cat("Calculating numbers of cell in Cluster", Cluster_label, "...\n")

  num <- data.frame(table(Template$cluster))
  num <- num[which(num$Var1 == Cluster_label), 2]

  return(num)
}

# 2.4: Cluster_sim() calculate the similarity of cells in cluster
#' Calculate the similarity of cells in cluster
#'
#' @param Template A table result contains cellIDs and cluster column.
#' @param SimCells  A matrix recorded similarity value between cell to cell.
#' @param Cluster_label A list of cells in specific cluster label.
#'
#' @return An integer represent the similarity in specific cluster.
#' @export
#'
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#'
#' @examples
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

# Function to calculate distance for a chunk of the matrix
#' Calculate the similarity value between cell to cell
#'
#' @param binsMatrix A matrix recorded each cell copy number in bin-level.
#'
#' @return A matrix recorded similarity value between cell to cell.
#' @export
#'
#' @examples
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
