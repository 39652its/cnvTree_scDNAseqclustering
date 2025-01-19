# 6.1: startTimed(), endTimed() for time calculating
startTimed <- function(...){
  x <- paste0(..., collapse = "")
  message(x, appendLF = FALSE)
  ptm <- proc.time()
  return(ptm)
}

endTimed <- function(ptm){
  time <- proc.time() - ptm
  message(" ", round(time[3], 2), "s")
}


# 6.2: QuickQC() based on Aneufinder clusterbyquality for a short cut
QuickQC <- function(Input, Spikiness, Bdistance){
  # Cluster by Quality
  message("Cluster by Quality ...")
  cl01_QCsb <- AneuFinder::clusterByQuality(Input, measures=c('spikiness','bhattacharyya', "entropy", "num.segments", "sos"))
  Cluster_parameters <- data.frame(cl01_QCsb$parameters)

  #存符合的組
  message("Selected cells ...")
  chosegroup <- subset(Cluster_parameters, spikiness<Spikiness & bhattacharyya>Bdistance)
  chosegroup <- as.numeric(rownames(chosegroup))

  selected <- lapply(chosegroup, function(index) cl01_QCsb$classification[[index]])
  selected <- unlist(selected)


  return(selected)
}


# 6.3: GenomeHeatmap() function for plotting CN pattern in each cells
GenomeHeatmap <- function(Input, cellID, pqArm_file){
  # import CN template and total cell copy number matrix
  CN_bins_template <- CN_template(input = Input, pqArm_file = pqArm_file)
  CN_chr_template <- CN_bins_template %>%
    group_by(chr) %>%
    summarise(length = max(end) - min(start)) %>%
    mutate(X_cum = c(0, cumsum(as.numeric(length[-nrow(.)]))))

  # Import the HC clustering cell order
  cellOrder <- clusterbyHMM(input = Input, selected = cellID)
  cellOrder <- cellOrder$IDorder %>%
    as.data.frame(.) %>%
    rownames_to_column(., var = "cellID")
  Input <- Input[cellOrder$cellID]

  # Import copy number data
  numofcell <- length(cellID)
  mat_input <- lapply(1:numofcell, function(i) segment_transform(Input[[i]], i, CN_chr_template))
  mat_input <- bind_rows(mat_input)
  mat_input$CN[which(mat_input$CN>5)] <- 5

  # plot chr site
  Chr_tmp <- CN_chr_template %>%
    mutate(chr = sub("^chr", "", chr),
           X_cum_end = length + X_cum,
           text_pos = (length + X_cum*2)/2,
           Y_end = numofcell-0.8)


  # Fixed color template: Copy number more than 5 as same color
  color_set <- c("#D0CECE", "#8165A3", "#9BBB59", "#C0504D", "#ED7D31", "#FFC000")
  del_neu_num <- which(as.numeric(levels(factor(mat_input$CN))) %in% c(0:5))
  color_tem <- c(color_set[del_neu_num])

  # figure legend
  lgd_labels <- sort(unique(mat_input$CN))
  lgd_labels[which(lgd_labels==5)] <- ">= 5"

  height = case_when(
    numofcell > 60 ~ numofcell,
    numofcell <= 60 & numofcell > 10 ~ 50,
    numofcell <= 10  ~ 30
  )

  PlotCN_heatmap <-
    ggplot(mat_input) +
    geom_segment(aes(x = X_cum_start, xend = X_cum_end, y = Y_cum_start, yend = Y_cum_end, color = factor(CN)), linewidth = 2, show.legend = TRUE) +
    scale_color_manual(name="Copy Number", labels=lgd_labels, values=color_tem)+
    geom_segment(data = Chr_tmp, aes(x = X_cum_end , y = -0.2, xend = X_cum_end, yend = Y_end)) +
    geom_text(data = Chr_tmp, aes(x = text_pos, y = -0.05 * height, label = chr), size = 10) +
    labs(x = "Chromosome", y = "Cell ID", color = "Copy Number") +
    xlim(c(0, max(Chr_tmp$X_cum_end)+1)) +
    ylim(c(-0.05 *height, numofcell)) +
    theme(plot.margin = margin(3, 2, 3, 1, "cm"),
          plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white")) +
    theme_void()


  return(PlotCN_heatmap)
}

# 6.3.1: Transform each cell CN-segment information
segment_transform <- function(data, index, CN_chr_template) {
  as.data.frame(data$bins) %>%
    mutate(segment = cumsum(copy.number != lag(copy.number, default = first(copy.number)))) %>%
    group_by(seqnames, segment, copy.number) %>%
    summarize(
      start = first(start),
      end = last(end),
      width = end - start + 1,
      strand = first(strand),
      .groups = "drop"
    ) %>%
    ungroup() %>%
    select(-c(segment, strand)) %>%
    rename("chr" = "seqnames", "CN" = "copy.number") %>%
    left_join(CN_chr_template, by = "chr") %>%
    mutate(
      X_cum_start = start + X_cum,
      X_cum_end = end + X_cum,
      Y_cum_start = index - 1,
      Y_cum_end = index - 1
    )
}

# 6.4: Totalcluster_pdf() function for creating pdf in total clusters by CN matrix
Totalcluster_pdf <- function(Input, Template, FILEname, FILEpath, pqArm_file, cellnum_name, cellcutoff, cluster_name){
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
    Fig_seq[[count]] <- GenomeHeatmap(Input = Input, cellID = SS, pqArm_file = pqArm_file)

  }
  Label <- paste0(Cluster_No, " (n = ", cellnum_list, ")")


  combine_plot <- ggarrange(plotlist = Fig_seq,
                            ncol= 1,
                            labels =Label,
                            common.legend = FALSE,
                            legend = "right",
                            hjust  = 0.8,
                            align = "v",
                            font.label = list(size = 45, face = "bold", color ="black"),
                            heights = c(height_ratio))+
    theme(plot.margin = margin(2,2,2,8, "cm"))



  ggsave(combine_plot,
         filename = paste0(FILEpath, FILEname),
         height = (200*sum(height_ratio) + 800*length(height_ratio)),
         width = 13000,
         units = "px",
         limitsize = FALSE)

  cat("Output ", FILEname, " is done. ","\n")
}


# 6.4.1: Totalcluster_Cluster_No() function for transfer the number of clusters in the data
Totalcluster_Cluster_No <- function(Template, cellnum_name, cellnum, cluster_name){
  Cluster_No <- Template %>%
    filter(get(cellnum_name) >= cellnum) %>%
    pull(get(cluster_name)) %>%
    unique(.) %>%
    sort(.)

  return(Cluster_No)
}

# 6.4.2: Totalcluster_SS() function for transfer each cluster of cells in the data
Totalcluster_SS <- function(Template, cluster_name, k){
  SS <- Template %>%
    filter(get(cluster_name) == k) %>%
    pull(cellID)

  return(SS)
}


# 6.5: writeOutput() function for outputting cluster results in .txt
writeOutput <- function(data, filename, path){
  FILEpath <- paste0(path, filename, ".txt")
  write.table(data, file = FILEpath, row.names = FALSE, col.names = TRUE)
}


# 6.6: scDNA.superimpose() function for superimpose between definedCNVs and DNA segment info
scDNA.superimpose <- function(Template, DefinedCNVs){
  Groups <- unique(Template$Subclone_CN$Subclone)
  Template$Subclone_CN$CNV_state <- Template$Subclone_CN$CN
  # CN only seperate in 3 types: del/neu/amp
  Template$Subclone_CN$CN = case_when(
    Template$Subclone_CN$CN <  2 ~ "del",
    Template$Subclone_CN$CN == 2 ~ "neu",
    Template$Subclone_CN$CN >  2 ~ "amp")

  superimpose <- NULL
  for(groups in 1:length(Groups)){
    intersection <- NULL
    # check defined CNVs in each group
    CNVs <- Template$Subclone_CN %>%
      filter(Subclone %in% c(Groups[groups]),
             chr %in% c(unique(DefinedCNVs$chr)),
             CN %in% c(DefinedCNVs$CN))

    intersection <- merge(DefinedCNVs, CNVs, by = c("chr", "CN"))
    intersection <- intersection %>%
      mutate(
        F_start = case_when(
          start < CNV_start ~ 0,
          start >= CNV_start & start <= CNV_end ~ 1,
          start > CNV_end ~ 2),
        F_end = case_when(
          end < CNV_start ~ 0,
          end >= CNV_start & end <= CNV_end ~ 1,
          end > CNV_end ~ 2),
        seg = paste0(F_start, F_end),
        final_start = case_when(
          seg %in% c("00", "22") ~ NA,
          seg %in% c("01", "02") ~ CNV_start,
          seg %in% c("11", "12") ~ start),
        final_end = case_when(
          seg %in% c("00", "22") ~ NA,
          seg %in% c("01", "11") ~ end,
          seg %in% c("02", "12") ~ CNV_end)) %>%
      filter(!seg %in% c("00", "22")) %>%
      mutate(cnv_range = final_end - final_start + 1) %>%
      group_by(CNV_region) %>%
      summarise(cnv_range = sum(cnv_range)) %>%
      as.data.frame(.)

    superimpose <- left_join(DefinedCNVs, intersection, by = "CNV_region") %>%
      mutate(Subclone = Groups[groups],
             CNV_range = CNV_end - CNV_start + 1,
             cnv_range = ifelse(is.na(cnv_range)==TRUE, 0, cnv_range),
             cnv_ratio = cnv_range / CNV_range,
             final_cnv = ifelse(cnv_ratio >= 0.5, 1, 0)) %>%
      rbind(., superimpose)
  }

  Template$superimpose <- superimpose


  return(Template)
}


# 6.7: scDNA.clustering() function for receiving DNA clustering output in superimpose range
scDNA.clustering <- function(Template){
  cnv_region <- unique(Template$superimpose$CNV_region)
  cnv_region <- paste0("CNV", cnv_region)
  Subclone_ss <- unique(Template$superimpose$Subclone)

  cnv_matrix <- NULL
  for(subclone in 1:length(Subclone_ss)){
    Cell_num <- Template$final_cluster_output %>%
      filter(Subclone %in% Subclone_ss[subclone]) %>%
      nrow(.)
    CNVs <- Template$superimpose %>%
      filter(Subclone %in% Subclone_ss[subclone]) %>%
      select(final_cnv) %>%
      rbind(Subclone_ss[subclone], Cell_num[1])
    if(subclone == 1){
      cnv_matrix <- CNVs
    } else {
      cnv_matrix <- cbind(cnv_matrix, CNVs)
    }
  }

  cnv_matrix <- as.data.frame(t(cnv_matrix)) %>%
    setNames(c(cnv_region, "DNA_cluster", "DNA_Cellnum")) %>%
    `rownames<-`(1:length(Subclone_ss))

  Template$DNA_cluster <- cnv_matrix


  return(Template)
}

# 6.8: cluster w/ or w/o CNV pattern plot
# Input = scRNA_output
scDNA_CNVpattern <- function(Input, DeterminedCNVs, FILEpath, FILEname){
  Data <- Input %>% select(!c(DNA_Cellnum, DNA_cluster))
  Data <- as.matrix(as.data.frame(Data))
  rownames(Data) <- paste0(Input$DNA_cluster, " (n=", Input$DNA_Cellnum, ")")
  colnames(Data) <- DeterminedCNVs %>%
    mutate(chr = sub("^chr", "", DeterminedCNVs$chr),
           CN_type = ifelse(CN=="amp", "+", "-"),
           chr_cytoband = paste0(CN_type, " ", chr, "(", first_band, "-", last_band, ")")) %>%
    pull(chr_cytoband)

  col_fun <- c("#EFEDF6", "#E6B745")

  # h <- min(length(Input$RNA_cluster)*75, )
  # w <- min((length(Input)-3)*150, 20000)

  png(filename = paste0(FILEpath, FILEname),
      width = (length(Input))*100,
      height = (length(Input$DNA_cluster))*150)

  Oncoscan <- ComplexHeatmap::Heatmap(Data,
                                      name = "CNV exist",
                                      cluster_rows = TRUE,
                                      cluster_columns = FALSE,
                                      show_column_dend = FALSE,
                                      show_row_dend = TRUE,

                                      show_heatmap_legend = TRUE,


                                      col = col_fun,
                                      rect_gp = gpar(col = "white", lwd = 2),
                                      #cell_fun = cell_fun,
                                      row_names_side = c("right"),
                                      column_names_side = c("bottom"),


                                      row_names_gp = gpar(fontsize = 40),
                                      show_row_names = TRUE,
                                      row_title  = "scDNA clusters",
                                      row_title_gp = gpar(fontsize = 50, fontface = "bold"),
                                      row_dend_width = unit(5, "cm"),  # 调整行树状图的宽度
                                      row_dend_gp = gpar(lwd = 2, col = "black"),

                                      column_names_gp = gpar(fontsize = 30),
                                      column_title = "High-confidence CNV regions",
                                      column_title_gp = gpar(fontsize = 50, fontface = "bold"),
                                      column_title_side = "top",
                                      column_names_rot = 90,
                                      column_names_centered = FALSE,
                                      border = "black",

                                      heatmap_legend_param = list(
                                        title_gp = gpar(fontsize = 40, fontface = "bold"),  # 調整圖例標題的字體大小
                                        labels_gp = gpar(fontsize = 35),  # 調整圖例標示的字體大小
                                        grid_width = unit(3, "cm"),
                                        grid_height = unit(3, "cm")
                                      )

  )

  # lgd = Legend(labels = c("0", "1"),
  #              legend_gp = gpar(fill = col_fun),
  #              title = "CNV exist",
  #              title_gp = gpar(fontsize = 25),
  #              labels_gp = gpar(fontsize = 18),
  #              grid_width = unit(1.75, "cm"),
  #              grid_height = unit(1.5, "cm"))
  #


  draw(Oncoscan, padding = unit(c(5, 2, 2, 2), "cm"))

  dev.off()

}
