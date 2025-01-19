# 7.0.1: Quick QC
# 1. Aneufider
QuickQC_Aneufider <- function(input, Spikiness, Bdistance){
  # Cluster by Quality
  message("Cluster by Quality ...")
  cl01_QCsb <- clusterByQuality(input, measures=c('spikiness','bhattacharyya', "entropy", "num.segments",
                                                  "sos"))
  Cluster_parameters <- data.frame(cl01_QCsb$parameters)

  #存符合的組
  selected.files <- NULL
  chosegroup <- NULL
  list(chosegroup)
  print("Selected files ...")
  for(i in 1:nrow(Cluster_parameters)){
    if ((Cluster_parameters$spikiness[i] < Spikiness ) && (Cluster_parameters$bhattacharyya[i] > Bdistance)){
      chosegroup <- c(chosegroup, i)
    }
  }
  selected.files <- unlist(cl01_QCsb$classification[chosegroup])

  return(selected.files)
}

# 7.1: Clustering() stands for clustering method
Clustering <- function(input, Cellnum=10, Sim=0.9){
  message("=== Step 01: Clustering ===")
  message("Started calculating number of cells and similarity")

  # cluster template
  Clust_cuttree <- CutTree_final(input = input, selected = names(input))

  # cell similarity
  seq_identity <- CN_seq(input = input, Template = names(input))
  Sim_cell2cell <- Cluster_SimTem(binsMatrix = seq_identity)

  Clustering_01 <- TRUE
  meet_criteria_label <- NULL
  meet_criteria_sim <- NULL
  count = 0
  while(Clustering_01 == TRUE){
    count = count + 1
    message(sprintf("Step 01: Clustering %d time%s\n", count, ifelse(count > 1, "s", "")))

    Cluster_label <- unique(Clust_cuttree$cluster[!Clust_cuttree$cluster %in% meet_criteria_label])
    cat("Cluster_label", Cluster_label, "\n")

    All_Cell_num <- NULL
    All_similarity <- NULL
    # Number of cells
    All_Cell_num <- sapply(Cluster_label, function(i) Cluster_num(Template = Clust_cuttree, Cluster_label = i))
    Clust_cuttree <- Clust_cuttree %>% filter(!(cluster %in% Cluster_label[All_Cell_num < Cellnum]))

    # Similarity
    Cluster_label <- unique(Clust_cuttree$cluster[!Clust_cuttree$cluster %in% meet_criteria_label])
    for (i in Cluster_label){
      Similarity <- Cluster_sim(SimCells = Sim_cell2cell, Template = Clust_cuttree, Cluster_label = i)
      if (Similarity < Sim){
        Clust_cuttree <- CutTree(input = input, Template = Clust_cuttree, Cluster_label = i)
        All_similarity <- c(All_similarity, Similarity)
      } else {
        All_similarity <- c(All_similarity, Similarity)
        meet_criteria_label <- c(meet_criteria_label, i)
        meet_criteria_sim <- c(meet_criteria_sim, Similarity)
      }
    }

    # cat("All_Cell_num", All_Cell_num, "\n")
    # cat("All_similarity", All_similarity, "\n")
    # cat("meet_criteria_label", meet_criteria_label, "\n")
    # cat("meet_criteria_sim", meet_criteria_sim, "\n")

    # 確定是否要再重複while loop
    if(length(All_Cell_num) == 0){
      if(length(meet_criteria_label) == 0){
        message("Please lower the similarity or # cell citeria in each cluster")
        Clustering_01 = FALSE
      } else {
        Clustering_01 = FALSE
      }
    } else if (min(All_Cell_num) > Cellnum && min(All_similarity) > Sim){
      Clustering_01 = FALSE
    }
  }



  return(Clust_cuttree)
}


# 7.2: pqArmClustering() stands for pqArm clustering method
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


# 7.3: Reclustering() stands for Reclustering method
Reclustering <- function(input, pqArm_output, pqArm_file, difratio_chr=0.3){
  message("=== Step 03: Reclustering ===")

  Recluster_Output <- NULL
  for(Label in unique(pqArm_output$cluster)){
    ptm <- startTimed("Reclustering for cluster ", Label, " ...")
    pqArm_merge_target <- NULL
    New_pqArm_clustering <- pqArm_recluster(pqArm_cluster = pqArm_output, Cluster = Label)

    # potential merge target must more than one (with less10 clusters( >=2 & <10), or more than one more10 clusters(>2))
    if(length(New_pqArm_clustering)>1){
      pqArm_dif <- pqArm_reclustering_dif(input = Sample, pqArm_recluster_sim = New_pqArm_clustering,
                                          pqArm_cluster = pqArm_output, Cluster = Label, pqArm_file = pqArm_file)
      pqArm_merge_target <- pqArm_reclusterBy_ratio_target(pqArm_cluster = pqArm_output, Cluster = Label, pqReclsut_sim = pqArm_dif, difratio_chr = difratio_chr)
    }

    # is.null(pqArm_merge_target)
    if(is.null(pqArm_merge_target) == FALSE){
      Recluster_CellID <- pqArm_recluster_result(pqArm_cluster = pqArm_output, Cluster = Label, pqReclsut_target = pqArm_merge_target)
      Recluster_Output <- rbind(Recluster_Output, Recluster_CellID)
    } else {
      Recluster_CellID <- pqArm_output %>%
        filter(cluster %in% Label) %>%
        mutate(Recluster_pattern = pqArm_pattern,
               Recluster_cellnum = pqArm_pattern_cellnum,
               Recluster_cluster = pqArm_cluster,
               PQreturn = NA,
               dif_num = NA,
               dif_ratio = NA)
      Recluster_Output <- rbind(Recluster_Output, Recluster_CellID)
    }
    endTimed(ptm)
  }

  return(Recluster_Output)
}


# 7.4: SubcloneClustering() stands for Subclone clsutering method
SubcloneClustering <- function(input, Reclustering_output, min_cell=5, overlap_region=10**7, dif_ratio=0.2){
  message("=== Step 04: Subclone Clustering ===")

  Final_output <- list()
  Subclone_cluster <- NULL
  Subclone_CNr <- NULL
  cell_clustering <- NULL

  # from overlap_region to number of bins in coverage
  binsize <- as.data.frame(Sample[[1]]$bins)$width[1]
  overlap_bp <- round(overlap_region/binsize, digits = 0)

  # Select enough cell number Reclusters, avoiding keep filtering
  R <- Reclustering_output %>% filter(Recluster_cellnum >= min_cell)
  for (Recluster_label in unique(R$Recluster_cluster)){
    ptm <- startTimed("Subclone clustering for Recluster ", Recluster_label, " ... \n")

    breakpoints <- collect_cluster_bp(input = input, Clustering_output = R, Recluster_label = Recluster_label)
    Cell_num <- R %>% filter(Recluster_cluster == Recluster_label) %>% pull(Recluster_cellnum)
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

  Final_output <- list(final_cluster_output = left_join(Reclustering_output, Subclone_cluster, by = "cellID"),
                       Subclone_CN = Subclone_CNr)


  return(Final_output)
}




# 7.5: scDNA_Output() stands for outputting scDNA clustering final outputs
scDNA_Output <- function(input, Summary, pqArm_file, output.dir=getwd(), consecutive_region=10**7, cellcutoff=5){
  message("=== Step 05: Output cnvTree results ===")


  # 1. cellID summary
  writeOutput(data = Summary$final_cluster_output, filename = "/cnvTree.scDNAseq_grouping", path = output.dir)

  # 2. Each region copy number to subclone
  writeOutput(data = Summary$Subclone_CN, filename = "/cnvTree.scDNAseq_SubcloneRegionCN", path = output.dir)

  # 3. All CNV region in the sample
  CNV_Data <- Total_cnvRegion(input = Sample, Template = Summary$Subclone_CN, pqArm_file = pqArm_file, consecutive_region = consecutive_region)
  CNV_Data <- cnvRegion.toPQarm(FILE = CNV_Data, pqArm_file = pqArm_file)
  writeOutput(data = CNV_Data, filename = "/cnvTree.scDNAseq_DefinedCNVregion", path = output.dir)

  # 4. DNA superimpose
  DNA_superimpose <- scDNA.superimpose(Template = Summary, DefinedCNVs = CNV_Data)
  DNA_superimpose <- scDNA.clustering(Template = DNA_superimpose)
  writeOutput(data = DNA_superimpose$DNA_cluster, filename = "/cnvTree.scDNAseq_DNAcluster", path = output.dir)

  # 5. Final cluster heatmap
  Totalcluster_pdf(Input = input, Template = Summary$final_cluster_output,
                   FILEname = "/cnvTree.scDNAseq_Grouping_fig.pdf", FILEpath = output.dir, pqArm_file = pqArm_file,
                   cellnum_name = "Subclone_cellnum", cellcutoff = cellcutoff, cluster_name = "Subclone")

  # 6. cluster heatmap with dendrogram
  scDNA_CNVpattern(Input=DNA_superimpose$DNA_cluster, DeterminedCNVs=CNV_Data, FILEpath=output.dir, FILEname="/cnvTree.scDNAseq_heatmap.png")
}
