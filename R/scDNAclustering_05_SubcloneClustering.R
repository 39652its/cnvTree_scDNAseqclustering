# 5.0: collect_cluster_bp() collect all break points in the cluster
collect_cluster_bp <- function(input, Clustering_output, Recluster_label){
  selected_files <- Clustering_output %>%
    filter(Recluster_cluster %in% Recluster_label) %>%
    pull(cellID)

  # Use lapply to gather breakpoints for all selected files at once
  breakpoints_list <- lapply(selected_files, function(i) {
    input[[i]]$breakpoints %>%
      as.data.frame() %>%
      mutate(cellID = i)
  })

  # Combine all the results using bind_rows, which is more efficient than rbind in a loop
  breakpoints <- bind_rows(breakpoints_list)

  return(breakpoints)
}

# 5.1: output_bp_covers() calcuate each bp with cover bps
output_bp_covers <- function(Template, binsize, overlap=10, overlap_times=Cell_num[1]*0.5){
  chr <- unique(Template$seqnames) %>%
    factor(levels = levels(factor(Template$seqnames))) %>%
    sort()

  # binsize <- input$width[1]

  # 將10個bp轉換為實際數字
  bp <- NULL
  bp_list <- vector("list", length(chr))

  bp_list <- lapply(chr, function(k) {
    # 將位點從10**6轉換成bins的格式
    chr_unique <- Template %>%
      filter(seqnames %in% k) %>%
      mutate(site = round(start / binsize, digits = 0)) %>%
      group_by(site) %>%
      summarise(times = n(), .groups = 'drop') %>%
      mutate(site = as.numeric(as.character(site)))

    # 將其中點的數值做延伸並合併
    site_vec <- chr_unique$site
    cover_matrix <- sapply(site_vec, function(x) {
      cover_label <- site_vec %in% ((x - overlap):(x + overlap))
      cover_nums <- sum(cover_label)
      cover_times <- sum(chr_unique$times[cover_label])
      return(c(cover_nums, cover_times))
    })

    chr_unique <- chr_unique %>%
      mutate(cover_nums = as.integer(cover_matrix[1, ]),
             cover_times = as.integer(cover_matrix[2, ])) %>%
      arrange(desc(cover_nums))

    # 選擇符合條件的位點
    bp_order <- chr_unique %>%
      filter(cover_times >= overlap_times) %>%
      arrange(desc(cover_nums)) %>%
      pull(site)

    if (length(bp_order) > 0) {
      # 初始化bp_group
      bp_group <- data.frame()
      count <- 1
      while (length(bp_order) > 0) {
        select_label <- which(bp_order %in% ((bp_order[1] - overlap):(bp_order[1] + overlap)))
        bp_group_1 <- data.frame(breakpoints = bp_order[select_label], events = count)
        bp_group <- bind_rows(bp_group, bp_group_1)
        bp_order <- bp_order[-select_label]
        count <- count + 1
      }

      # 將結果存入列表
      bp_group <- bp_group %>%
        mutate(chr = k)
      bp_site_select <- merge(chr_unique, bp_group, by.x = "site", by.y = "breakpoints")
      return(bp_site_select)
    }
  })

  bp <- bind_rows(bp_list) %>%
    select(chr, site, times, events, cover_nums, cover_times) %>%
    mutate(binsize = binsize)

  cat("All the breakpoints in same Recluster group of cells ... \n")

  return(bp)
}

# 5.2: bp_events() use the distribution of breakpoints to define the events
bp_events <- function(input, Template, binsize){
  # set chr levels
  vec <- unique(Template$chr)
  nums <- as.numeric(gsub("chr", "", vec)[grepl("\\d", vec)])
  nums <- paste0("chr", nums[order(nums)])
  # Levels <- c(nums, vec[!grepl("\\d", vec)])
  Template$chr <- factor(Template$chr, levels = nums)  # 更加安全的做法是直接轉換因子

  # binsize <- Template$binsize[1]
  event_region <- NULL

  event_region_list <- lapply(nums, function(k) {
    event <- unique(Template$events[which(Template$chr == k)])

    # 使用 lapply 來處理每個 event
    event_region_1 <- lapply(event, function(i) {
      Template %>%
        filter(chr == k, events == i) %>%
        summarise(
          event = i,
          cross_bp = n(),
          min_site = min(site, na.rm = TRUE),
          max_site = max(site, na.rm = TRUE),
          chr = k
        )
    })

    # 使用 bind_rows 合併 event
    bind_rows(event_region_1)
  })

  # 合併所有的 chr 結果
  event_region <- bind_rows(event_region_list) %>%
    arrange(chr, min_site) %>%
    mutate(event = row_number())  # 使用 row_number() 確保 event 重新編號

  desired_order <- c("chr", "event", "cross_bp", "min_site", "max_site")

  event_region <- event_region[ ,desired_order] %>%
    mutate(chr = factor(chr, levels = nums)) %>%
    arrange(chr)

  event_region <- event_region.bin(input = input, Template = event_region)
  event_region$binsize <- binsize


  return(event_region)
}

# 5.4.1: event_region.bin() change event_region to bins level(sites)
event_region.bin <- function(input, Template){
  bins_num <- input[[1]]$bins@seqnames %>%
    table(.) %>%
    data.frame(.) %>%
    setNames(c("chr", "Freq"))

  bins_num$CDF_start <- sapply(1:nrow(bins_num), function(x){
    sum(bins_num$Freq[1:x-1])
  })
  bins_num$CDF_start <- bins_num$CDF_start+1
  bins_num$CDF_end <- sapply(1:nrow(bins_num), function(x){
    sum(bins_num$Freq[1:x])
  })

  Template <- left_join(bins_num, Template, by = "chr") %>%
    mutate(chr = factor(chr, levels = levels(bins_num$chr)),
           bins_level_min = (min_site + CDF_start - 1),
           bins_level_max = (max_site + CDF_start - 1),
           event = ifelse(is.na(event), 0, event))




  return(Template)
}

# 5.4: bp_region() from defined-event get the region sites
bp_region <- function(event, binsize){
  # binsize = event$binsize[1]
  region <- data.frame()
  for (i in levels(event$chr)){
    event_R <- event %>%
      filter(chr %in% c(i)) %>%
      mutate(region = event,
             bp_start = min_site*binsize,
             bp_end = max_site*binsize-1) %>%
      select(chr, bp_start, bp_end, region, Freq, CDF_start, CDF_end, bins_level_min, bins_level_max)

    start <- NULL
    end <- NULL
    event_binstart <- NULL
    event_binend <- NULL
    if(event_R$region[1] == 0){  # some chrmosome with no breakpoint
      event_R <- event_R %>%
        mutate(bp_start = 1,
               bp_end = event_R$Freq[1]*binsize,
               region_ratio = 1) %>%
        select(chr, bp_start, bp_end, region, Freq, region_ratio, CDF_start, CDF_end) %>%
        setNames(c("chr", "start", "end", "region", "region_size", "region_ratio", "event_binstart", "event_binend"))
    } else {
      for (j in 1:(nrow(event_R)+1)){
        if(j == 1){
          start <- c(start, 1)
          end <- c(end, event_R$bp_start[j]-1)
          event_binstart <- c(event_binstart, event_R$CDF_start[j])
          event_binend <- c(event_binend, event_R$bins_level_min[j])
        } else if (j == (nrow(event_R)+1)){
          start <- c(start, event_R$bp_end[j-1]+1)
          end <- c(end, event_R$Freq[j-1]*binsize)
          event_binstart <- c(event_binstart, event_R$bins_level_max[j-1])
          event_binend <- c(event_binend, event_R$CDF_end[j-1])
        } else {
          start <- c(start, event_R$bp_end[j-1]+1)
          end <- c(end, event_R$bp_start[j]-1)
          event_binstart <- c(event_binstart, event_R$bins_level_max[j-1])
          event_binend <- c(event_binend, event_R$bins_level_min[j])
        }
      }
      event_R <- event_R %>%
        add_row() %>%
        mutate(chr = lag(chr, default = first(as.character(chr))),
               CDF_start = lag(CDF_start , default = first(CDF_start)),
               CDF_end = lag(CDF_end , default = first(CDF_end)),
               region = 1:n(),
               start = start,
               end = end,
               event_binstart = event_binstart,
               event_binend = event_binend,
               region_size = (event_binend-event_binstart+1),
               region_ratio = (region_size/(CDF_end-CDF_start+1))) %>%
        select(chr, start, end, region, region_size, region_ratio, event_binstart, event_binend)
    }

    region <- event_R %>%
      rbind(region)

  }

  region <- region %>%
    mutate(chr = factor(chr, levels = levels(event$chr)),
           region_size = round(region_size, 0)) %>%
    arrange(chr)


  return(region)
}

# 5.5: Region_CN() smooth CN based on defined regions
Region_CN <- function(input, Reclustering_output, Recluster_label, events){

  selected_files <- Reclustering_output %>%
    filter(Recluster_cluster %in% Recluster_label) %>%
    pull(cellID)

  CN_matrix <- CN_seq(input = input, Template = selected_files)

  Smooth_CN <- data.frame()
  for(i in 1:nrow(events)){
    # cat("Dealing with chr", events$chr[i], " Region", events$region[i], "copy number......\n")
    R_binsCN <- CN_matrix[events$event_binstart[i]:events$event_binend[i],selected_files]
    R_CN <- sapply(1:ncol(R_binsCN), function(a){
      freq <- table(R_binsCN[ ,a]) %>%
        as.data.frame() %>%
        arrange(desc(Freq)) %>%
        pull(Var1) %>%
        as.character() %>%
        as.integer()
      first_element <- freq[1]
    })

    Smooth_CN <- Smooth_CN %>%
      rbind(R_CN)

  }
  Smooth_CN <- Smooth_CN %>%
    setNames(selected_files)

  return(Smooth_CN)

}

# 5.6: Subclone_clustering() check cell to cell whether with >?% of different chromosome than clustering
Subclone_clustering <- function(CN_incells_input, event_region, dif_ratio, Subclone_num){
  cell_code <- c(colnames(CN_incells_input))

  difChr_num <- NULL
  ptm <- startTimed("Calculating cell to cell different ratio ......") # Time code
  # Convert 'event_region' to data.table
  event_region_dt <- as.data.table(event_region)

  # Preallocate matrix for results
  difChr_num <- matrix(0, nrow = length(cell_code), ncol = length(cell_code))

  for (i in 1:(length(cell_code) - 1)) {
    for (j in (i + 1):length(cell_code)) {
      # Find the differing regions
      dif_region <- which(CN_incells_input[, i] != CN_incells_input[, j])

      # Aggregate differences using data.table for speed
      dif_chr_ratio <- event_region_dt[dif_region, .(total_region_ratio = sum(region_ratio)), by = chr][
        total_region_ratio > dif_ratio, .N]

      # Assign the result to both (i, j) and (j, i) due to symmetry
      difChr_num[i, j] <- dif_chr_ratio
      difChr_num[j, i] <- dif_chr_ratio
    }
  }
  endTimed(ptm)


  difChr_num <- as.data.frame(difChr_num) %>%
    setNames(cell_code)
  rownames(difChr_num) <- cell_code

  # cell to cell : Matrix about number of chromosomes
  Check_num <- sapply(1:ncol(difChr_num), function(x){
    str <- length(which(difChr_num[ , x] == 0))
  })

  # a <- difChr_num
  Subclone <- NULL
  count = Subclone_num
  Subclone <- Check_num %>%
    as.data.frame(.) %>%
    mutate(cellID = cell_code,
           subclone = NA) %>%
    setNames(c("Subclone_cellnum", "cellID", "Subclone")) %>%
    arrange(desc(Subclone_cellnum)) %>%
    select(c("cellID", "Subclone"))


  while(any(is.na(Subclone$Subclone))){
    count = count + 1
    ss <- min(which(is.na(Subclone$Subclone) == TRUE))
    ss <- Subclone$cellID[ss]
    selected <- which(difChr_num[ , ss] == 0)
    selected <- rownames(difChr_num)[selected]
    difChr_num <- difChr_num[!rownames(difChr_num) %in% selected, ]

    Subclone$Subclone[which(Subclone$cellID %in% selected)] <- count
  }

  Cellnum <- as.data.frame(table(Subclone$Subclone)) %>%
    setNames(c("Subclone", "Subclone_cellnum"))

  Subclone <- merge(Subclone, Cellnum, by = "Subclone")


  return(Subclone)
}

### Subclone cnv outputs
# 5.7: Subclone_CNregion() output each region copy number in each subclone, as the basement of CNV template
Subclone_CNregion <- function(sep_region, CN_region, each_subclone, min_cell, output=c("SubcloneCNVRegion", "SubcloneRegionCN")){
  s <- each_subclone %>%
    filter(Subclone_cellnum >= min_cell)
  Subclone_CN <- NULL
  for(Label in unique(s$Subclone)){
    ss <- s %>%
      filter(Subclone %in% Label) %>%
      pull(cellID)
    s_CN <- CN_region[ ,ss]
    R_CN <- sapply(1:nrow(s_CN), function(a){
      freq <- as.numeric(s_CN[a, ]) %>%
        table() %>%
        as.data.frame() %>%
        arrange(desc(Freq)) %>%
        setNames(c("CN", "Freq")) %>%
        pull(CN) %>%
        as.character() %>%
        as.integer()
      first_element <- freq[1]
    })

    Sub_CN <- sep_region %>%
      mutate(Subclone = Label,
             CN = R_CN) %>%
      select(c(chr, start, end, region, Subclone, CN))

    Subclone_CN <- Subclone_CN %>%
      rbind(Sub_CN)
  }

  # decide the output
  if(output == "SubcloneCNVRegion"){
    Subclone_CN <- Subclone_CN %>%
      filter(!CN %in% 2)
  } else if (output == "SubcloneRegionCN"){
    Subclone_CN <- Subclone_CN
  } else {
    message("ERROR: Not found the output")
  }

  return(Subclone_CN)
}

# 5.8: Total_cnvRegion() output total cnv regions across subclones
Total_cnvRegion <- function(input, Template, pqArm_file, consecutive_region){
  CN_tem <- data.frame(seqnames(input[[1]]$bins), ranges(input[[1]]$bins))
  CN_tem <- CN_tem %>%
    setNames(c("chr", "start", "end", "width"))

  # consecutive_bins
  consecutive_bins = round(consecutive_region/CN_tem$width[1], digits = 0)

  Final_CNVr <- list()
  Final_CNV <- NULL
  Final_CNVr$del <- Total_cnvRegion.DelAmp(Template = Template, CN_tem = CN_tem, method = c("Del"))
  Final_CNVr$amp <- Total_cnvRegion.DelAmp(Template = Template, CN_tem = CN_tem, method = c("Amp"))


  # Masked centromere region
  # pqArm_range <- pqArm_file.remake(FILE = pqArm_file)
  Masked <- pqArm_file.cen(FILE = pqArm_file)
  for(CN_type in names(Final_CNVr)){
    Final_CNVr[[CN_type]] <- Final_CNVr[[CN_type]] %>%
      merge(., Masked) %>%
      filter((start >= MaskEnd | end <= MaskStart),
             n_subclone != 0)

    Final_CNVr[[CN_type]] <- Final_CNVr[[CN_type]] %>%
      mutate(gap = cumsum(c(0, diff(rows) != 1)),
             chr_gapno = paste0(chr, "_",gap)) %>%
      group_by(chr_gapno) %>%
      filter(n() > consecutive_bins) %>%
      summarise(CNV_start = min(start),
                CNV_end = max(end)) %>%
      separate(chr_gapno, into = c("chr", "CNV_region"), sep = "_") %>%
      as.data.frame()

    Final_CNVr[[CN_type]] <- Final_CNVr[[CN_type]] %>% mutate(CN = CN_type)
  }
  Final_CNV <- Final_CNVr$del %>%
    rbind(Final_CNVr$amp) %>%
    mutate(chr = factor(chr, levels = levels(CN_tem$chr))) %>%
    arrange(chr) %>%
    mutate(CNV_region = 1:n()) %>%
    select(c("chr", "CNV_region", "CN", "CNV_start", "CNV_end"))


  return(Final_CNV)
}

# 5.8.1: Total_cnvRegion() calculate the cnv happends in Deletion and Amplification in whole chromosome
Total_cnvRegion.DelAmp <- function(Template, CN_tem, method = c("Del", "Amp")){
  if (method == "Del"){
    s_chr <- Template %>% filter(CN < 2)
  } else if (method == "Amp"){
    s_chr <- Template %>% filter(CN > 2)
  }

  s_chr <- s_chr %>%
    arrange(chr) %>%
    select(chr) %>%
    unique(.) %>%
    pull(.)

  # check whole chromosome region each site with occurs how many subclone CNVs
  Final_CNVr <- NULL
  for(k in s_chr){
    ss <- CN_tem %>%
      filter(chr %in% k) %>%
      mutate(rows = 1:nrow(.),
             n_subclone = 0)

    if (method == "Del"){
      S <- Template %>%
        filter(CN < 2, chr%in%k) %>%
        arrange(chr, CN, start)
    } else if (method == "Amp"){
      S <- Template %>%
        filter(CN > 2, chr%in%k) %>%
        arrange(chr, CN, start)
    }
    for (i in 1:nrow(S)){
      Selected <- ss %>%
        filter(start >= S$start[i], end <= S$end[i]) %>%
        pull(rows)
      ss$n_subclone[Selected] <- ss$n_subclone[Selected]+1
    }
    Final_CNVr <- rbind(Final_CNVr, ss)
  }

  return(Final_CNVr)
}

# 5.8.2: cnvRegion.toPQarm() input cnvRegions in defined format then add p, q arm information automatically
cnvRegion.toPQarm <- function(FILE, pqArm_file){
  pqArm_range <- read.table(gzfile(pqArm_file),sep="\t", col.names = c("chr", "start", "end", "name","gieStain")) %>%
    filter(chr %in% FILE$chr)

  Intersect <- merge(FILE, pqArm_range, by = "chr") %>%
    mutate(final_start = case_when(CNV_start<start ~ 1,
                                   CNV_start>=start & CNV_start<=end ~ 2,
                                   CNV_start>end ~ 3),
           final_end = case_when(CNV_end<start ~ 1,
                                 CNV_end>=start & CNV_end<=end ~ 2,
                                 CNV_end>end ~ 3),
           pattern = paste0(final_start, final_end),
           overlap = ifelse(pattern %in% c(11, 33), 0, 1)) %>%
    filter(overlap == 1) %>%
    arrange(chr, CNV_region  , start) %>%
    group_by(CNV_region)  %>%
    summarise(first_band = first(name), last_band = last(name))

  Final_output <- merge(FILE, Intersect, by = "CNV_region")


  return(Final_output)
}
