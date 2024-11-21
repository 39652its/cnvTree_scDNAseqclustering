##### 01: Input
# 1.0: changeFormat() read .rds to change as GRanges format
#' Change .rds copy number data as GRanges format
#'
#' @param file A table with column names: "cellID", "seqnames", "start", "end", and "copy.number" in .rds format.
#' @param cores A integer. parallel execution in number of cores, default=1.
#'
#' @return a list with each cell in GRanges object.
#' @export
#'
#' @examples
#' file_path <- system.file("extdata", "example_data.rds", package = "cnvTree")
#' Example_data <- changeFormat(file = file_path, core = 4)
#'
changeFormat <- function(file, cores) {
  ptm <- startTimed("Read RDS file...")

  Bin_CN <- readRDS(file)
  data.table::setDT(Bin_CN)
  endTimed(ptm)

  # Set up parallel environment
  future::plan(future::multisession, workers = cores)
  ptm <- startTimed("Make GRanges object ...")


  # Split data by "cellID"
  Bin_CN_list <- split(Bin_CN, by = "cellID", keep.by = FALSE)
  template_data <- Bin_CN_list[[1]][, setdiff(names(Bin_CN_list[[1]]), "copy.number"), with = FALSE]
  template_gr <- GenomicRanges::makeGRangesFromDataFrame(template_data, keep.extra.columns = TRUE)

  NewFormat <- list()
  NewFormat <- lapply(seq_along(Bin_CN_list), function(i) {
    # Extract current Bin_CN_list element
    Bins <- template_gr
    cell_data <- Bin_CN_list[[i]]

    # Copy the template GRanges object and add "copy.number"
    Bins <- template_gr
    S4Vectors::mcols(Bins)$copy.number <- cell_data$copy.number

    # Compute breakpoints
    shifted_cn <- data.table::shift(cell_data$copy.number, type = "lag", fill = NA)
    breakpoints <- cell_data[copy.number != shifted_cn, .(seqnames, start, end, copy.number)]

    # Convert breakpoints to GRanges if not empty
    Breakpoints <- if (nrow(breakpoints) == 0) {
      NULL
    } else {
      GenomicRanges::GRanges(
        seqnames = breakpoints$seqnames,
        ranges = IRanges::IRanges(start = breakpoints$start, end = breakpoints$end),
        copy.number = breakpoints$copy.number
      )
    }

    # Return a list with the new format
    list(ID = names(Bin_CN_list)[i], bins = Bins, breakpoints = Breakpoints)
  })

  names(NewFormat) <- names(Bin_CN_list)
  endTimed(ptm)

  return(NewFormat)
}


# 1.1: get the copy number matrix
#' Gather copy number matrix in selected cells
#'
#' @param input a list with each cell in GRanges object.
#' @param Template a list of selected cells.
#'
#' @return a integer matrix, which column names are selected cellIDs and row names are genome ranges separated in fixed bins.
#' @export
#'
#' @examples
#' file_path <- system.file("extdata", "example_data.rds", package = "cnvTree")
#' Example_data <- changeFormat(file = file_path, core = 4)
#' cellIDs <- names(Example_data)[1:10]
#' CNmatrix <- CN_seq(input = Example_data, Template = cellIDs)
#'
CN_seq <- function(input, Template){

  c <- list()
  c <- lapply(Template, function(i) {
    as.integer(input[[i]][["bins"]]$copy.number)}
  )

  c <- as.data.frame(c)
  colnames(c) <- Template

  return(c)
}

