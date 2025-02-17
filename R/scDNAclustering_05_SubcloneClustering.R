# 5.0: collect_cluster_bp() collect all break points in the cluster
#' Collect breakpoints from a cluster of cells
#'
#' This function extracts all breakpoints from a specified cluster of cells
#' based on the results of the re-clustering step.
#'
#' @param input A named list where each element is a `GRanges` object representing a single cell.
#' @param Clustering_output A data frame recording the pqArm clustering,
#'   and re-clustering results for each cell. This table tracks the clustering history at each step.
#' @param Recluster_label An integer specifying the cluster from the re-clustering step
#'   for which breakpoints should be extracted.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @return A data frame containing breakpoint sites with the following columns:
#'   - `seqnames`: Chromosome name (chr1, chr2, ...).
#'   - `start`: Start position of the breakpoint.
#'   - `end`: End position of the breakpoint.
#'   - `width`: Width of the fixed-bin size in the genomic regions.
#'   - `strand`: Strand information (`+` or `-`).
#'   - `copy.number`: Left copy number value at the breakpoint.
#'   - `cellID`: Identifier of the corresponding cell.
#'
#' @keywords internal
#'
collect_cluster_bp <- function(input, Clustering_output, Recluster_label){
  selected_files <- Clustering_output %>%
    dplyr::filter(.data$Recluster_cluster %in% Recluster_label) %>%
    dplyr::pull(.data$cellID)

  # Use lapply to gather breakpoints for all selected files at once
  breakpoints_list <- lapply(selected_files, function(i) {
    input[[i]]$breakpoints %>%
      as.data.frame() %>%
      dplyr::mutate(cellID = i)
  })

  # Combine all the results using bind_rows, which is more efficient than rbind in a loop
  breakpoints <- dplyr::bind_rows(breakpoints_list)

  return(breakpoints)
}

# 5.1: output_bp_covers() calcuate each bp with cover bps
#' Identify significant breakpoints based on coverage in a specified range
#'
#' This function calculates how many breakpoints can be covered within a specified range
#' for each breakpoint, treating it as the midpoint. It also filters significant breakpoints
#' based on a defined threshold.
#'
#' @param Template A data frame containing breakpoint sites with the following columns:
#'   - `seqnames`: Chromosome name (chr1, chr2, ...).
#'   - `start`: Start position of the breakpoint.
#'   - `end`: End position of the breakpoint.
#'   - `width`: Width of the fixed-bin size in the genomic regions.
#'   - `strand`: Strand information (`+` or `-`).
#'   - `copy.number`: Left copy number value at the breakpoint.
#'   - `cellID`: Identifier of the corresponding cell.
#' @param binsize An integer specifying the fixed-bin size, referring to the copy number
#'   variation calling result where the genomic region is split into bins.
#' @param overlap An integer defining the copy number frequently changed regions.
#' @param overlap_times A numeric value used as a filtering criterion for significant breakpoints in event region calculation.
#'   The default value is set as half the number of cells in the cluster.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @return A data frame recording significant breakpoints with the following columns:
#'   - `chr`: Chromosome name (chr1, chr2, ...).
#'   - `site`: Breakpoint position.
#'   - `times`: The frequency of breakpoint occurrence.
#'   - `events`: Number of distinct breakpoint events.
#'   - `cover_nums`: Number of breakpoints covered within the specified range.
#'   - `cover_times`: Number of times a breakpoint is covered.
#'   - `binsize`: Fixed-bin size used in the calculation.
#'
#' @keywords internal
#'
output_bp_covers <- function(Template, binsize, overlap, overlap_times){
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
      dplyr::filter(.data$seqnames %in% k) %>%
      dplyr::mutate(site = round(.data$start / binsize, digits = 0)) %>%
      dplyr::group_by(.data$site) %>%
      dplyr::summarise(times = dplyr::n(), .groups = 'drop') %>%
      dplyr::mutate(site = as.numeric(as.character(.data$site)))

    # 將其中點的數值做延伸並合併
    site_vec <- chr_unique$site
    cover_matrix <- sapply(site_vec, function(x) {
      cover_label <- site_vec %in% ((x - overlap):(x + overlap))
      cover_nums <- sum(cover_label)
      cover_times <- sum(chr_unique$times[cover_label])
      return(c(cover_nums, cover_times))
    })

    chr_unique <- chr_unique %>%
      dplyr::mutate(cover_nums = as.integer(cover_matrix[1, ]),
                    cover_times = as.integer(cover_matrix[2, ])) %>%
      dplyr::arrange(dplyr::desc(.data$cover_nums))

    # 選擇符合條件的位點
    bp_order <- chr_unique %>%
      dplyr::filter(.data$cover_times >= overlap_times) %>%
      dplyr::arrange(dplyr::desc(.data$cover_nums)) %>%
      dplyr::pull(.data$site)

    if (length(bp_order) > 0) {
      # 初始化bp_group
      bp_group <- data.frame()
      count <- 1
      while (length(bp_order) > 0) {
        select_label <- which(bp_order %in% ((bp_order[1] - overlap):(bp_order[1] + overlap)))
        bp_group_1 <- data.frame(breakpoints = bp_order[select_label], events = count)
        bp_group <- dplyr::bind_rows(bp_group, bp_group_1)
        bp_order <- bp_order[-select_label]
        count <- count + 1
      }

      # 將結果存入列表
      bp_group <- bp_group %>%
        dplyr::mutate(chr = k)
      bp_site_select <- merge(chr_unique, bp_group, by.x = "site", by.y = "breakpoints")
      return(bp_site_select)
    }
  })

  bp <- dplyr::bind_rows(bp_list) %>%
    dplyr::select(.data$chr, .data$site, .data$times, .data$events, .data$cover_nums, .data$cover_times) %>%
    dplyr::mutate(binsize = binsize)

  cat("All the breakpoints in same Recluster group of cells ... \n")

  return(bp)
}

# 5.2: bp_events() use the distribution of breakpoints to define the events
#' Define copy number variation hotspots based on significant breakpoints
#'
#' This function identifies genomic regions where copy number variations (CNVs)
#' frequently occur, based on a list of significant breakpoints.
#'
#' @param input A named list where each element is a `GRanges` object representing a single cell.
#' @param Template A data frame recording significant breakpoints with the following columns:
#'   - `chr`: Chromosome name (chr1, chr2, ...).
#'   - `site`: Breakpoint position.
#'   - `times`: The frequency of breakpoint occurrence.
#'   - `events`: Number of distinct breakpoint events.
#'   - `cover_nums`: Number of breakpoints covered within the specified range.
#'   - `cover_times`: Number of times a breakpoint is covered.
#'   - `binsize`: Fixed-bin size used in the calculation.
#' @param binsize An integer specifying the fixed-bin size, referring to the copy number
#'   variation calling result where the genomic region is split into bins.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @return A data frame listing CNV events for each chromosome, identifying regions where copy number variations frequently occur.
#'
#' @keywords internal
#'
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
        dplyr::filter(.data$chr == k, .data$events == i) %>%
        dplyr::summarise(
          event = i,
          cross_bp = dplyr::n(),
          min_site = min(.data$site, na.rm = TRUE),
          max_site = max(.data$site, na.rm = TRUE),
          chr = k
        )
    })

    # 使用 bind_rows 合併 event
    dplyr::bind_rows(event_region_1)
  })

  # 合併所有的 chr 結果
  event_region <- dplyr::bind_rows(event_region_list) %>%
    dplyr::arrange(.data$chr, .data$min_site) %>%
    dplyr::mutate(event = dplyr::row_number())  # 使用 row_number() 確保 event 重新編號

  desired_order <- c("chr", "event", "cross_bp", "min_site", "max_site")

  event_region <- event_region[ ,desired_order] %>%
    dplyr::mutate(chr = factor(.data$chr, levels = nums)) %>%
    dplyr::arrange(.data$chr)

  event_region <- event_region.bin(input = input, Template = event_region)
  event_region$binsize <- binsize


  return(event_region)
}

# 5.4.1: event_region.bin() change event_region to bins level(sites)
#' Convert event region data to continuous bin-level sites
#'
#' This function transforms event region data from discrete chromosome sites
#' into continuous bin-level sites for further analysis.
#'
#' @param input A named list where each element is a `GRanges` object representing a single cell.
#' @param Template A data frame recording copy number variation (CNV) events for each chromosome.
#'   This table will be updated with bin-level site information.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @return A data frame `Template` with additional columns:
#'   - `bins_level_min`: The minimum bin-level position for each event region.
#'   - `bins_level_max`: The maximum bin-level position for each event region.
#'   These columns enable further bin-level calculations.
#'
#' @keywords internal
#'
event_region.bin <- function(input, Template){
  bins_num <- input[[1]]$bins@seqnames %>%
    table() %>%
    data.frame() %>%
    stats::setNames(c("chr", "Freq"))

  bins_num$CDF_start <- sapply(1:nrow(bins_num), function(x){
    sum(bins_num$Freq[1:x-1])
  })
  bins_num$CDF_start <- bins_num$CDF_start+1
  bins_num$CDF_end <- sapply(1:nrow(bins_num), function(x){
    sum(bins_num$Freq[1:x])
  })

  Template <- dplyr::left_join(bins_num, Template, by = "chr") %>%
    dplyr::mutate(chr = factor(.data$chr, levels = levels(bins_num$chr)),
                  bins_level_min = (.data$min_site + .data$CDF_start - 1),
                  bins_level_max = (.data$max_site + .data$CDF_start - 1),
                  event = ifelse(is.na(.data$event), 0, .data$event))




  return(Template)
}

# 5.4: bp_region() from defined-event get the region sites
#' Define segments in continuous bin-level sites based on event regions
#'
#' This function segments the genome into continuous bin-level regions based on an event region template.
#'
#' @param event A data frame recording event regions in continuous bin-level sites,
#'   where each row represents an individual event region.
#' @param binsize An integer specifying the fixed-bin size, referring to the copy number
#'   variation calling result where the genomic region is split into bins.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @return A data frame representing chromosome segments at a continuous bin-level resolution,
#'   divided according to the event region template.
#'
#' @keywords internal
#'
bp_region <- function(event, binsize){
  # binsize = event$binsize[1]
  region <- data.frame()
  for (i in levels(event$chr)){
    event_R <- event %>%
      dplyr::filter(.data$chr %in% c(i)) %>%
      dplyr::mutate(region = .data$event,
                    bp_start = .data$min_site*binsize,
                    bp_end = .data$max_site*binsize) %>%
      dplyr::select(.data$chr, .data$bp_start, .data$bp_end, .data$region, .data$Freq,
                    .data$CDF_start, .data$CDF_end, .data$bins_level_min, .data$bins_level_max)

    start <- NULL
    end <- NULL
    event_binstart <- NULL
    event_binend <- NULL
    if(event_R$region[1] == 0){  # some chrmosome with no breakpoint
      event_R <- event_R %>%
        dplyr::mutate(bp_start = 1,
                      bp_end = event_R$Freq[1]*binsize,
                      region_ratio = 1) %>%
        dplyr::select(.data$chr, .data$bp_start, .data$bp_end, .data$region, .data$Freq, .data$region_ratio, .data$CDF_start, .data$CDF_end) %>%
        stats::setNames(c("chr", "start", "end", "region", "region_size", "region_ratio", "event_binstart", "event_binend"))
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
        tibble::add_row() %>%
        dplyr::mutate(chr = dplyr::lag(.data$chr, default = gdata::first(as.character(.data$chr))),
                      CDF_start = dplyr::lag(.data$CDF_start , default = gdata::first(.data$CDF_start)),
                      CDF_end = dplyr::lag(.data$CDF_end , default = gdata::first(.data$CDF_end)),
                      region = 1:dplyr::n(),
                      start = start,
                      end = end,
                      event_binstart = event_binstart,
                      event_binend = event_binend,
                      region_size = (.data$event_binend-.data$event_binstart+1),
                      region_ratio = (.data$region_size/(.data$CDF_end-.data$CDF_start+1))) %>%
        dplyr::select(.data$chr, .data$start, .data$end, .data$region, .data$region_size, .data$region_ratio, .data$event_binstart, .data$event_binend)
    }

    region <- event_R %>%
      rbind(region)

  }

  region <- region %>%
    dplyr::mutate(chr = factor(.data$chr, levels = levels(event$chr)),
                  region_size = round(.data$region_size, 0)) %>%
    dplyr::arrange(.data$chr)


  return(region)
}

# 5.5: Region_CN() smooth CN based on defined regions
#' Calculate copy number for chromosome segments
#'
#' This function calculates the copy number in each chromosome segment based on
#' a continuous bin-level segment template. It uses re-clustering results and event
#' region templates to determine copy number variations (CNVs) for each cell.
#'
#' @param input A named list where each element is a `GRanges` object representing a single cell.
#' @param Reclustering_output A data frame recording the pqArm clustering,
#' and re-clustering step results for each cell, including the clustering history.
#' @param Recluster_label An integer specifying the cluster from the re-clustering step
#' for which copy number calculation is performed.
#' @param events A data frame containing continuous bin-level chromosome segments, divided by event region templates.
#'
#' @return A numeric data frame where each column corresponds to a cell and each row
#' maps to a chromosome segment from the input event region template.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @keywords internal
#'
Region_CN <- function(input, Reclustering_output, Recluster_label, events){

  selected_files <- Reclustering_output %>%
    dplyr::filter(.data$Recluster_cluster %in% Recluster_label) %>%
    dplyr::pull(.data$cellID)

  CN_matrix <- CN_seq(input = input, Template = selected_files)

  Smooth_CN <- data.frame()
  for(i in 1:nrow(events)){
    # cat("Dealing with chr", events$chr[i], " Region", events$region[i], "copy number......\n")
    R_binsCN <- CN_matrix[events$event_binstart[i]:events$event_binend[i],selected_files]
    R_CN <- sapply(1:ncol(R_binsCN), function(a){
      freq <- table(R_binsCN[ ,a]) %>%
        as.data.frame() %>%
        dplyr::arrange(dplyr::desc(.data$Freq)) %>%
        dplyr::pull(.data$Var1) %>%
        as.character() %>%
        as.integer()
      first_element <- freq[1]
    })

    Smooth_CN <- Smooth_CN %>%
      rbind(R_CN)

  }
  Smooth_CN <- Smooth_CN %>%
    stats::setNames(selected_files)

  return(Smooth_CN)

}

# 5.6: Subclone_clustering() check cell to cell whether with >?% of different chromosome than clustering
#' Perform subclone clustering based on copy number patterns
#'
#' This function performs subclone clustering by analyzing copy number patterns in predefined chromosome segment templates.
#' It groups cells into subclones by comparing copy number variations (CNVs) while tolerating a specified level
#' of difference between cells.
#'
#' @param CN_incells_input A numeric data frame where each column corresponds to a cell,
#' and each row maps to a chromosome segment from the segment template.
#' @param event_region A data frame containing continuous bin-level chromosome segments,
#' divided by event region templates.
#' @param dif_ratio A numeric value specifying the tolerance threshold for differences
#' in copy number patterns between cells during clustering.
#' @param Subclone_num An integer used as the index for the subclone.
#'
#' @return A data frame recording the final result of subclone clustering, containing:
#'   - `Subclone`: The subclone identifier.
#'   - `cellID`: The unique identifier for each cell.
#'   - `Subclone_cellnum`: The number of cells in each subclone.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @keywords internal
#'
Subclone_clustering <- function(CN_incells_input, event_region, dif_ratio, Subclone_num){
  cell_code <- c(colnames(CN_incells_input))

  difChr_num <- NULL
  ptm <- startTimed("Calculating cell to cell different ratio ......") # Time code
  # Convert 'event_region' to data.table
  event_region_dt <- data.table::as.data.table(event_region)

  # Preallocate matrix for results
  difChr_num <- matrix(0, nrow = length(cell_code), ncol = length(cell_code))

  for (i in 1:(length(cell_code) - 1)) {
    for (j in (i + 1):length(cell_code)) {
      # Find the differing regions
      dif_region <- which(CN_incells_input[, i] != CN_incells_input[, j])

      # Aggregate differences using data.table for speed
      dif_event_region <- event_region_dt[dif_region, ]
      dif_chr_ratio <- dif_event_region %>%
        dplyr::group_by(.data$chr) %>%
        dplyr::summarise(total_region_ratio = sum(.data$region_ratio), .groups = "drop") %>%
        dplyr::filter(.data$total_region_ratio > dif_ratio) %>%
        dplyr::summarise(count = dplyr::n()) %>%
        dplyr::pull(count)

      if (length(dif_chr_ratio) == 0) {
        dif_chr_ratio <- 0
      }


      # Assign the result to both (i, j) and (j, i) due to symmetry
      difChr_num[i, j] <- dif_chr_ratio
      difChr_num[j, i] <- dif_chr_ratio
    }
  }
  endTimed(ptm)


  difChr_num <- as.data.frame(difChr_num) %>%
    stats::setNames(cell_code)
  rownames(difChr_num) <- cell_code

  # cell to cell : Matrix about number of chromosomes
  Check_num <- sapply(1:ncol(difChr_num), function(x){
    str <- length(which(difChr_num[ , x] == 0))
  })

  # a <- difChr_num
  Subclone <- NULL
  count = Subclone_num
  Subclone <- Check_num %>%
    as.data.frame() %>%
    dplyr::mutate(cellID = cell_code,
                  subclone = NA) %>%
    stats::setNames(c("Subclone_cellnum", "cellID", "Subclone")) %>%
    dplyr::arrange(dplyr::desc(.data$Subclone_cellnum)) %>%
    dplyr::select(c("cellID", "Subclone"))


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
    stats::setNames(c("Subclone", "Subclone_cellnum"))

  Subclone <- merge(Subclone, Cellnum, by = "Subclone")


  return(Subclone)
}

### Subclone cnv outputs
# 5.7: Subclone_CNregion() output each region copy number in each subclone, as the basement of CNV template
#' Output copy number by unique segment template for each subclone
#'
#' This function generates copy number (CN) outputs for each subclone using a unique chromosome segment template.
#' It filters clusters based on a minimum cell count and outputs CN values for each region.
#'
#' @param sep_region A data frame recording chromosome segments, divided by event region template.
#' @param CN_region A numeric data frame where each column corresponds to a cell,
#' and each row maps to a chromosome segment from the segment template.
#' @param each_subclone A data frame recording the clustering history, including pqArm clustering,
#' re-clustering, and subclone clustering results for each cell.
#' @param min_cell An integer specifying the minimum cell count required for a cluster to be included in the output.
#' @param output A character string specifying the column name in `each_subclone` to be used as output information.
#'
#' @return A data frame recording each subclone's unique chromosome segment template and its corresponding copy number (CN).
#'  The output includes the following columns:
#'   - `chr`: Chromosome name (chr1, chr2, ...).
#'   - `start`: Start position of the chromosome segment.
#'   - `end`: End position of the chromosome segment.
#'   - `region`: Unique region identifier.
#'   - `Subclone`: Subclone identifier.
#'   - `CN`: Copy number for the corresponding segment.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' \dontrun{
#' subclone_cn_output <- Subclone_CNregion(
#'   sep_region = event_template,
#'   CN_region = copy_number_matrix,
#'   each_subclone = subclone_results,
#'   min_cell = 10,
#'   output = "Subclone"
#' )
#' head(subclone_cn_output)
#' }
Subclone_CNregion <- function(sep_region, CN_region, each_subclone, min_cell, output=c("SubcloneCNVRegion", "SubcloneRegionCN")){
  s <- each_subclone %>%
    dplyr::filter(.data$Subclone_cellnum >= min_cell)
  Subclone_CN <- NULL
  for(Label in unique(s$Subclone)){
    ss <- s %>%
      dplyr::filter(.data$Subclone %in% Label) %>%
      dplyr::pull(.data$cellID)
    s_CN <- CN_region[ ,ss]
    R_CN <- sapply(1:nrow(s_CN), function(a){
      freq <- as.numeric(s_CN[a, ]) %>%
        table() %>%
        as.data.frame() %>%
        dplyr::arrange(dplyr::desc(.data$Freq)) %>%
        stats::setNames(c("CN", "Freq")) %>%
        dplyr::pull(.data$CN) %>%
        as.character() %>%
        as.integer()
      first_element <- freq[1]
    })

    Sub_CN <- sep_region %>%
      dplyr::mutate(Subclone = Label,
                    CN = R_CN) %>%
      dplyr::select(c(.data$chr, .data$start, .data$end, .data$region, .data$Subclone, .data$CN))

    Subclone_CN <- Subclone_CN %>%
      rbind(Sub_CN)
  }

  # decide the output
  if(output == "SubcloneCNVRegion"){
    Subclone_CN <- Subclone_CN %>%
      dplyr::filter(!.data$CN %in% 2)
  } else if (output == "SubcloneRegionCN"){
    Subclone_CN <- Subclone_CN
  } else {
    message("ERROR: Not found the output")
  }

  return(Subclone_CN)
}

# 5.8: Total_cnvRegion() output total cnv regions across subclones
#' Output defined CNV regions across clusters
#'
#' This function identifies and outputs copy number variation (CNV) regions across subclones,
#' using a predefined segment template and cytoband information.
#' The CNV regions are classified into two types: amplifications ("amp") and deletions ("del").
#'
#' @param input A named list where each element is a `GRanges` object representing a single cell.
#' @param Template A data frame recording each subclone's unique chromosome segment template and its corresponding copy number (CN).
#'  The output includes the following columns:
#'   - `chr`: Chromosome name (chr1, chr2, ...).
#'   - `start`: Start position of the chromosome segment.
#'   - `end`: End position of the chromosome segment.
#'   - `region`: Unique region identifier.
#'   - `Subclone`: Subclone identifier.
#'   - `CN`: Copy number for the corresponding segment.
#' @param pqArm_file A table for cytoband information seen on Giemsa-stained chromosomes.
#' It should include the following columns:
#'   - `chrom`: Reference sequence chromosome or scaffold.
#'   - `chromStart`: Start position in genoSeq.
#'   - `chromEnd`: End position in genoSeq.
#'   - `name`: Name of cytogenetic band.
#'   - `gieStain`: Giemsa stain results.
#' @param consecutive_region A numeric value specifying the minimum length criteria for filtering CNV regions, default=10^7 bp.
#'
#'
#' @return A data frame recording CNV regions across subclones, containing the following columns:
#'   - `chr`: Chromosome name (chr1, chr2, ...).
#'   - `CNV_region`: The index of CNV regions.
#'   - `CN`: Copy number state, categorized as either "amp" (Amplification) or "del" (Deletion).
#'   - `CNV_start`: Start position of the CNV region.
#'   - `CNV_end`: End position of the CNV region.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' \dontrun{
#' cnv_regions <- Total_cnvRegion(
#'   input = single_cell_data,
#'   Template = subclone_template,
#'   pqArm_file = cytoband_data,
#'   consecutive_region = 3
#' )
#' head(cnv_regions)
#' }
#'
Total_cnvRegion <- function(input, Template, pqArm_file, consecutive_region){
  CN_tem <- data.frame(GenomicRanges::seqnames(input[[1]]$bins), IRanges::ranges(input[[1]]$bins))
  CN_tem <- CN_tem %>%
    stats::setNames(c("chr", "start", "end", "width"))

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
      merge(Masked) %>%
      dplyr::filter((.data$start >= .data$MaskEnd | .data$end <= .data$MaskStart),
                    .data$n_subclone != 0)

    Final_CNVr[[CN_type]] <- Final_CNVr[[CN_type]] %>%
      dplyr::mutate(gap = cumsum(c(0, diff(.data$rows) != 1)),
                    chr_gapno = paste0(.data$chr, "_", .data$gap)) %>%
      dplyr::group_by(.data$chr_gapno) %>%
      dplyr::filter(dplyr::n() > consecutive_bins) %>%
      dplyr::summarise(CNV_start = min(.data$start),
                       CNV_end = max(.data$end)) %>%
      tidyr::separate(.data$chr_gapno, into = c("chr", "CNV_region"), sep = "_") %>%
      as.data.frame()

    Final_CNVr[[CN_type]] <- Final_CNVr[[CN_type]] %>% dplyr::mutate(CN = CN_type)
  }
  Final_CNV <- Final_CNVr$del %>%
    rbind(Final_CNVr$amp) %>%
    dplyr::mutate(chr = factor(.data$chr, levels = levels(CN_tem$chr))) %>%
    dplyr::arrange(.data$chr) %>%
    dplyr::mutate(CNV_region = seq_len(dplyr::n())) %>%
    dplyr::select(c("chr", "CNV_region", "CN", "CNV_start", "CNV_end"))


  return(Final_CNV)
}

# 5.8.1: Total_cnvRegion() calculate the cnv happends in Deletion and Amplification in whole chromosome
#' Classify copy number variations as deletion or amplification
#'
#' This function identifies and classifies copy number variations (CNVs) into two types:
#' Deletion ("Del") or Amplification ("Amp"). It uses a fixed-bin size genome template
#' and processes each CNV type separately to generate bin-level summaries.
#'
#' @param Template A data frame recording each subclone's unique chromosome segment template and its corresponding copy number (CN).
#'  The output includes the following columns:
#'   - `chr`: Chromosome name (chr1, chr2, ...).
#'   - `start`: Start position of the chromosome segment.
#'   - `end`: End position of the chromosome segment.
#'   - `region`: Unique region identifier.
#'   - `Subclone`: Subclone identifier.
#'   - `CN`: Copy number for the corresponding segment.
#' @param CN_tem A data frame with columns "chr", "start", "end", and "width", representing the genome divided into fixed-bin size segments.
#' @param method A character string specifying the CNV type to calculate. Must be one of type `Del`, or `Amp`, corresponding to Deletion or Amplification.
#'
#' @return A list with two data frames: `del` and `amp`. Each table contains the following columns:
#'   - `chr`: Chromosome name (chr1, chr2, ...).
#'   - `start`: Start position of the bin in the CNV region.
#'   - `end`: End position of the bin in the CNV region.
#'   - `width`: Width of the bin.
#'   - `rows`: Row indices corresponding to the bins.
#'   - `n_subclone`: Number of subclones with CNV in the bin.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @keywords internal
#'
Total_cnvRegion.DelAmp <- function(Template, CN_tem, method = c("Del", "Amp")){
  if (method == "Del"){
    s_chr <- Template %>% dplyr::filter(.data$CN < 2)
  } else if (method == "Amp"){
    s_chr <- Template %>% dplyr::filter(.data$CN > 2)
  }

  s_chr <- s_chr %>%
    dplyr::arrange(.data$chr) %>%
    dplyr::select(.data$chr) %>%
    base::unique(.) %>%
    dplyr::pull(.)

  # check whole chromosome region each site with occurs how many subclone CNVs
  Final_CNVr <- NULL
  for(k in s_chr){
    ss <- CN_tem %>%
      dplyr::filter(.data$chr %in% k) %>%
      dplyr::mutate(rows = seq_len(nrow(.)),
                    n_subclone = 0)

    if (method == "Del"){
      S <- Template %>%
        dplyr::filter(.data$CN < 2,
                      .data$chr%in%k) %>%
        dplyr::arrange(.data$chr, .data$CN, .data$start)
    } else if (method == "Amp"){
      S <- Template %>%
        dplyr::filter(.data$CN > 2, .data$chr%in%k) %>%
        dplyr::arrange(.data$chr, .data$CN, .data$start)
    }
    for (i in 1:nrow(S)){
      Selected <- ss %>%
        dplyr::filter(.data$start >= S$start[i], .data$end <= S$end[i]) %>%
        dplyr::pull(.data$rows)
      ss$n_subclone[Selected] <- ss$n_subclone[Selected]+1
    }
    Final_CNVr <- rbind(Final_CNVr, ss)
  }

  return(Final_CNVr)
}

# 5.8.2: cnvRegion.toPQarm() input cnvRegions in defined format then add p, q arm information automatically
#' Map copy number variation regions to cytoband sites
#'
#' This function maps copy number variation (CNV) regions to cytoband sites,
#' based on Giemsa-stained chromosome information. It annotates CNV regions
#' with the corresponding cytoband labels for the first and last affected bands.
#'
#' @param FILE Either a character string specifying the file path for output,
#' or a connection open for writing. An empty string (\code{""}) indicates output to the console.
#' @param pqArm_file A table for cytoband information seen on Giemsa-stained chromosomes.
#' It should include the following columns:
#'   - `chrom`: Reference sequence chromosome or scaffold.
#'   - `chromStart`: Start position in genoSeq.
#'   - `chromEnd`: End position in genoSeq.
#'   - `name`: Name of cytogenetic band.
#'   - `gieStain`: Giemsa stain results.
#'
#'
#' @return A data frame mapping CNV regions to cytoband sites, containing the following columns:
#'   - `chr`: Chromosome name (chr1, chr2, ...).
#'   - `CNV_region`: The index of CNV regions.
#'   - `CN`: Copy number state, categorized as either "amp" (Amplification) or "del" (Deletion).
#'   - `CNV_start`: Start position of the CNV region.
#'   - `CNV_end`: End position of the CNV region.
#'   - `first_band`: Cytoband label of the first affected band.
#'   - `last_band`: Cytoband label of the last affected band.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @keywords internal
#'
cnvRegion.toPQarm <- function(FILE, pqArm_file){
  pqArm_range <- utils::read.table(gzfile(pqArm_file),sep="\t", col.names = c("chr", "start", "end", "name","gieStain")) %>%
    dplyr::filter(.data$chr %in% FILE$chr)

  Intersect <- merge(FILE, pqArm_range, by = "chr") %>%
    dplyr::mutate(final_start = dplyr::case_when(.data$CNV_start<.data$start ~ 1,
                                                 .data$CNV_start>=.data$start & .data$CNV_start<=.data$end ~ 2,
                                                 .data$CNV_start>.data$end ~ 3),
                  final_end = dplyr::case_when(.data$CNV_end<.data$start ~ 1,
                                               .data$CNV_end>=.data$start & .data$CNV_end<=.data$end ~ 2,
                                               .data$CNV_end>.data$end ~ 3),
                  pattern = paste0(.data$final_start, .data$final_end),
                  overlap = ifelse(.data$pattern %in% c(11, 33), 0, 1)) %>%
    dplyr::filter(.data$overlap == 1) %>%
    dplyr::arrange(.data$chr, .data$CNV_region, .data$start) %>%
    dplyr::group_by(.data$CNV_region)  %>%
    dplyr::summarise(first_band = gdata::first(.data$name),
                     last_band = gdata::last(.data$name))

  Final_output <- merge(FILE, Intersect, by = "CNV_region")


  return(Final_output)
}
