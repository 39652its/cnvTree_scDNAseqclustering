##### 02: Clustering

# 2.0: clusterbyHMM() based on AneuFinder::clusterHMMs(), calcuate distance than hierrachical clustering
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
CutTree_final <- function(input, selected){
  # 分群的原始檔，後面要用他作為基底
  message("Clustering and Data processing ...")

  clust <- clusterbyHMM(input = input, selected = selected)
  Clust_cuttree <- data.frame(cluster = cutree(clust[["hclust"]], k = 2),
                              cell = names(clust$IDorder))

  Clust_cuttree <- Clust_cuttree %>%
    mutate(cluster = as.numeric(cluster)) %>%
    arrange(desc(cluster))


  return(Clust_cuttree)
}

# 2.2: CutTree() Cut the tree repeatly
CutTree <- function(input, Template, Cluster_label){
  # selected.files建立
  # Cluster_label: Clust_cuttree$cluster中的分群數字
  selected.files <- NULL
  selected.files <- subset(Template, cluster%in%c(Cluster_label))$cell

  message("Divided the cluster in k=2 ")
  clust <- clusterbyHMM(input = input, selected = selected.files, exclude.regions = NULL)
  Clust_2 <- data.frame(NewCluster = cutree(clust[["hclust"]], k = 2),
                        cell = names(clust$IDorder))


  Template <- merge(Template, Clust_2, by = "cell", all = TRUE)
  Template$NewCluster <- ifelse(is.na(Template$NewCluster)==T, 0, Template$NewCluster)

  Clust_2 <- which(is.na(Template$NewCluster) == F)

  count = max(Template$cluster, na.rm = TRUE)

  Template$cluster <- case_when(
    Template$NewCluster == 1 ~ count + 1,
    Template$NewCluster == 2 ~ count + 2,
    TRUE ~ Template$cluster
  )

  Template <- subset(Template, select = -NewCluster)

  return(Template)
}

# 2.3: Cluster_num() calculate the number of cells in cluster
Cluster_num <- function(Template, Cluster_label){
  # cat("Calculating numbers of cell in Cluster", Cluster_label, "...\n")

  num <- data.frame(table(Template$cluster))
  num <- num[which(num$Var1 == Cluster_label), 2]

  return(num)
}

# 2.4: Cluster_sim() calculate the similairty of cells in cluster
Cluster_sim <- function(Template, SimCells, Cluster_label){
  # cat("Calculating cell similarity in Cluster ",  Cluster_label, " ...\n")

  selected <- Template %>% filter(cluster %in% Cluster_label) %>% pull(cell)
  selected <- which(SimCells$IDorder %in% selected)

  Similarity <- SimCells$similarity[selected, selected]
  Similarity[is.na(Similarity)] <- 0
  Similarity <- mean(Similarity)


  return(Similarity)
}

# Function to calculate distance for a chunk of the matrix
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
