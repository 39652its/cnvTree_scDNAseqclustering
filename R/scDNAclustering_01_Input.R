##### 01: Input
# 1.0: changeFormat() read .rds to change as GRanges format
changeFormat <- function(file, cores) {
  ptm <- startTimed("Read RDS file...")

  Bin_CN <- readRDS(file)
  setDT(Bin_CN)
  endTimed(ptm)

  # 设置并行环境
  plan(multisession, workers = cores)
  ptm <- startTimed("Make GRanges object ...")


  # 将数据按照 cellID 分组，形成列表
  Bin_CN_list <- split(Bin_CN, by = "cellID", keep.by = FALSE)
  Tem <- Bin_CN_list[[1]] %>% select(!copy.number)
  Tem <- makeGRangesFromDataFrame(Tem, keep.extra.columns = TRUE)

  NewFormat <- list()
  NewFormat <- lapply(seq_along(Bin_CN_list), function(i) {
    # Extract current Bin_CN_list element
    Bins <- Tem
    mcols(Bins)$copy.number <- Bin_CN_list[[i]]$copy.number

    # Compute breakpoints
    breakpoints <- Bin_CN_list[[i]][, .SD[copy.number != data.table::shift(copy.number, n=1)], by = seqnames]

    # Create GRanges object for breakpoints
    Breakpoints <- if (nrow(breakpoints) == 0) {
      NULL
    } else {
      GRanges(seqnames = breakpoints$seqnames,
              ranges = IRanges(start = breakpoints$start, end = breakpoints$end),
              copy.number = breakpoints$copy.number)
    }

    # Return a list with the new format
    list(ID = names(Bin_CN_list)[i], bins = Bins, breakpoints = Breakpoints)
  })

  names(NewFormat) <- names(Bin_CN_list)
  endTimed(ptm)

  return(NewFormat)
}
