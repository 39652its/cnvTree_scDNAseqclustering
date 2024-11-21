##### 04: Reclustering

# 4.1: pqArm_recluster() calculate the similarity between small clusters and >2 cells clusters
pqArm_recluster <- function(pqArm_cluster, Cluster){
  # pqArm_cluster <- read.xlsx(xlsxFile = FILEpath) ##這裡需要檢查資料的function
  pqArm_cluster <- pqArm_cluster %>%
    filter(cluster%in%Cluster)

  # unique pattern output
  Pattern_more10 <- pqArm_cluster %>%
    filter(pqArm_pattern_cellnum >= 2) %>%
    pull(pqArm_pattern) %>%
    unique()

  Pattern_less10 <- pqArm_cluster %>%
    filter(pqArm_pattern_cellnum < 10, pqArm_pattern_cellnum >= 2) %>%
    pull(pqArm_pattern) %>%
    unique()

  # check any Pattern_more10 or Pattern_less10 is NULL
  if(length(Pattern_more10)==0 | length(Pattern_less10)==0){
    return(0)
  } else{
    Pattern_more10 <- pqArm_cluster.pattern(Pattern_more10)
    Pattern_less10 <- pqArm_cluster.pattern(pattern = Pattern_less10)
  }


  # New_cluster similarity calculation
  New_cluster <- NULL
  N_cluster <- NULL
  for (i in 1:ncol(Pattern_more10)){
    for (j in 1:ncol(Pattern_less10)){
      N_cluster <- c(N_cluster, euclidean(Pattern_more10[ ,i], Pattern_less10[ ,j]))
    }
    New_cluster <- rbind(New_cluster, N_cluster)
    N_cluster <- NULL
  }

  New_cluster <- New_cluster %>%
    as.data.frame() %>%
    `row.names<-`(colnames(Pattern_more10)) %>%
    `colnames<-`(colnames(Pattern_less10))

  return(New_cluster)

}


# 4.1.1: pqArm_cluster.pattern() convert each cluster pattern from vector to sequence
pqArm_cluster.pattern <- function(pattern){
  P <- pattern %>%
    strsplit(. , split = "_") %>%
    data.frame() %>%
    setNames(c(pattern)) %>%
    mutate_if(is.character, as.numeric)

  return(P)
}


# 4.1.2: euclidean() euclidian distance calculation
euclidean <- function(a, b){
  sqrt(sum((a - b)^2))
}


# 4.2: pqArm_reclustering_dif() output the different ratio in different pqArm at bins-level between two clusters
pqArm_reclustering_dif <- function(input, pqArm_recluster_sim, pqArm_cluster, Cluster, pqArm_file){
  pqArm_sim <- apply(pqArm_recluster_sim, 2, function(x) min(x[x!=0]) )   #2: column is the more10 cluster pattern


  pqArm_sim <- as.data.frame(pqArm_sim)

  CN_bins_template <- CN_template(input = input, pqArm_file = pqArm_file)
  chr_pq <- paste0(rep(levels(CN_bins_template$chr), each = 2), rep(c("p", "q"), 11))

  # 將細胞數<10的cluster 與細胞數>2的cluster 相比，找到相似度最高的>2 Cluster 考慮合併
  new_pqArm_cluster <- NULL
  for(i in 1:nrow(pqArm_sim)){
    p <- rownames(pqArm_sim)[i]
    more10_Row <- which (pqArm_recluster_sim[,p] == pqArm_sim[i, 1])
    less10 <- rep(p, times = length(more10_Row))
    new_pqArm_cluster <- cbind(less10, rownames(pqArm_recluster_sim)[more10_Row]) %>%
      rbind(new_pqArm_cluster)
  }
  new_pqArm_cluster <- new_pqArm_cluster %>%
    as.data.frame() %>%
    setNames(c("less10", "more10"))

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
    filter(cluster == Cluster)
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
pqArm_return.PQ <- function(pattern, PQarm){
  pqArm_list <- PQarm

  Pattern_unlist <- pqArm_cluster.pattern(pattern = pattern) %>%
    setNames(c("less10", "more10"))

  pqArm_select <- which(Pattern_unlist$less10 != Pattern_unlist$more10) %>%
    pqArm_list[.]


  return(pqArm_select)
}


# 4.2.2: pqArm_return.Bins() select different pqArm to output the region at bin-level
pqArm_return.Bins <- function(Pattern, which_Arm, Tem, CN_matrix){
  ID <- Tem %>%
    filter(pqArm_pattern %in% c(Pattern)) %>%
    pull(cellID)

  CNmatrix <- CN_matrix %>%
    filter(Chr_arm %in% which_Arm) %>%
    subset(., select = c(ID))
  CNmatrix <- pqArm_DelNeuAmp(matrix = CNmatrix) #只看Del/Neu/Amp

  CN_bins <- sapply(1:nrow(CNmatrix), function(x) {
    freq <- table(as.integer(CNmatrix[x, ]))
    sorted_freq <- sort(freq, decreasing = TRUE)
    first_element <- as.integer(names(sorted_freq)[1])
    return(first_element)
  })


  return(CN_bins)
}


# 4.3: pqArm_reclusterBy_ratio_target() filter ratio and merge the clusters if criteria meets
pqArm_reclusterBy_ratio_target <- function(pqArm_cluster, Cluster, pqReclsut_sim, difratio_chr){
  # pqArm_cluster <- read.xlsx(xlsxFile = FILEpath)
  pqArm_cluster <- pqArm_cluster %>%
    filter(cluster %in% Cluster)

  Chioce <- pqReclsut_sim %>%
    group_by(less10) %>%
    mutate(less10_times = n(),
           merge_pattern = paste0(less10, "_", more10))
  cellnum <- pqArm_cluster %>%
    select(pqArm_pattern, pqArm_pattern_cellnum) %>%
    filter(pqArm_pattern %in% c(Chioce$more10)) %>%
    distinct(pqArm_pattern, .keep_all = TRUE) %>%
    setNames(c("more10", "more10_cellnum"))
  Chioce <- merge(Chioce, cellnum, by = "more10")


  # 多個region 不同的要都符合才能留下
  for (pattern in unique(Chioce$merge_pattern)){
    Selected <- Chioce %>%
      filter(merge_pattern %in% c(pattern),
             dif_ratio > difratio_chr)

    if(nrow(Selected)>0){
      Chioce <- Chioce %>%
        filter(!merge_pattern %in% c(pattern))
    } else {
      Chioce <- Chioce
    }
  }

  # 一種less10 最終只能配對到一個more10
  Chioce_Result <- NULL
  #pattern = unique(Chioce$less10)[3]
  for (pattern in unique(Chioce$less10)){
    Selected <- Chioce %>%
      select(less10, more10, PQreturn, dif_num, dif_ratio, more10_cellnum) %>%
      filter(less10 %in% c(pattern))
    if(length(unique(Selected$more10))>1){
      Selected <- Selected %>%
        filter(dif_ratio == min(Selected$dif_ratio))
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
    # 處理merge到的對象是cellnum<10的情況，因為這些群會同時出現在左邊與右邊
    more10_pattern_list <- Chioce_Result %>%
      filter(more10_cellnum<10) %>%
      arrange(desc(more10_cellnum)) %>%
      pull(more10)
    new_Chioceless10 <- NULL
    while (length(more10_pattern_list)>0){
      pattern =  more10_pattern_list[1]
      check_pattern <- Chioce_Result %>%
        filter(less10 == pattern | more10 == pattern)
      cellnum_list <- check_pattern %>%
        filter(more10_cellnum > 10)

      if (nrow(cellnum_list) > 0){
        new_pattern <- cellnum_list %>% pull(more10)
        new_Chioceless10 <- Chioce_Result %>%
          filter(less10 == pattern | more10 == pattern) %>%
          mutate(more10 = new_pattern) %>%
          rbind(new_Chioceless10)
      } else {
        new_Chioceless10 <- Chioce_Result %>%
          filter(less10 == pattern | more10 == pattern) %>%
          mutate(more10 = pattern) %>%
          rbind(new_Chioceless10)
      }
      rm_pattern <- unlist(unique(new_Chioceless10$less10))
      intersection <- intersect(more10_pattern_list, rm_pattern)
      more10_pattern_list <- setdiff(more10_pattern_list, intersection)

    }

    # 先將重新編輯過分群目的的細胞群移除，再併入最終的分群結果中
    Chioce_Result <- Chioce_Result %>%
      filter(!less10 %in% c(new_Chioceless10$less10)) %>%
      rbind(new_Chioceless10)

    return(Chioce_Result)
  }
}


# 4.4: pqArm_recluster_result() reset the content in pqArmCluster_CellID.xlsx and merge final result
pqArm_recluster_result <- function(pqArm_cluster, Cluster, pqReclsut_target){
  # pqArm_cluster <- read.xlsx(xlsxFile = FILEpath)
  pqArm_cluster = pqArm_cluster %>%
    filter(cluster %in% Cluster)

  pqArm_cluster <- left_join(pqArm_cluster, pqReclsut_target, by = c("pqArm_pattern" = "less10")) %>%
    select(!more10_cellnum) %>%
    dplyr::rename("Recluster_pattern" = "more10" )
  pqArm_cluster$Recluster_pattern <- ifelse(is.na(pqArm_cluster$Recluster_pattern) == TRUE, pqArm_cluster$pqArm_pattern, pqArm_cluster$Recluster_pattern)
  Recluster_summary <- pqArm_reclustering_summary(Data = pqArm_cluster$Recluster_pattern)
  pqArm_cluster <- left_join(pqArm_cluster, Recluster_summary , by = "Recluster_pattern") %>%
    arrange(desc(Recluster_cellnum)) %>%
    select(cellID, cluster, pqArm_pattern, pqArm_pattern_cellnum, pqArm_cluster, Recluster_pattern, Recluster_cellnum, Recluster_cluster, PQreturn, dif_num, dif_ratio)

  return(pqArm_cluster)
}


# 4.4.1: pqArm_reclustering_summary() create pqArm clustering final results
pqArm_reclustering_summary <- function(Data){
  cluster_table <- table(Data) %>%
    as.data.frame() %>%
    arrange(desc(Freq)) %>%
    mutate(Recluster_cluster = c(1:nrow(.))) %>%
    # filter(Freq > 1) %>%
    setNames(c("Recluster_pattern", "Recluster_cellnum", "Recluster_cluster"))

  return(cluster_table)
}

