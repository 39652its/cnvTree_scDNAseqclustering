##### 01: Input
# 1.0: changeFormat() read .rds or .txt to change as GRanges format
#' Convert copy number data to GRanges Format
#'
#' This function reads a `.rds` file or `.txt` file containing copy number variation (CNV) data
#' and converts it into a list of `GRanges` objects, where each element corresponds to a single cell.
#'
#' @param file A `.rds` file or `.txt` file containing a data frame with the following required columns:
#'   - `cellID`: Unique identifier for each cell.
#'   - `seqnames`: Chromosome or sequence name.
#'   - `start`: Start position of the segment.
#'   - `end`: End position of the segment.
#'   - `copy.number`: Copy number value for the segment.
#' @param cores An integer specifying the number of CPU cores to use for parallel processing, default=1.
#'
#' @return A named list where each element represents a cell, containing its corresponding genomic segments as a `GRanges` object.
#' @export
#'
#' @examples
#'
#' file_path <- system.file("extdata", "example_data.rds", package = "cnvTree")
#' Example_data <- changeFormat(file = file_path, core = 4)
#'
#'
changeFormat <- function(file, cores=1) {
  ptm <- startTimed("Read input file...")
  if(endsWith(file, ".rds") == TRUE){
    Bin_CN <- readRDS(file)
  } else if(endsWith(file, ".txt") == TRUE){
    Bin_CN <- utils::read.delim2(file, sep = " ")
  } else {
    "Please input either .rds or .txt format as input."
    endTimed(ptm)
    stop(changeFormat)
  }

  data.table::setDT(Bin_CN)
  endTimed(ptm)

  # Set up parallel environment
  future::plan(future::multisession, workers = cores)
  ptm <- startTimed("Make GRanges object ...")


  # Split data by "cellID"
  Bin_CN_list <- split(Bin_CN, by = "cellID", keep.by = FALSE)
  template_data <- Bin_CN_list[[1]][, setdiff(names(Bin_CN_list[[1]]), "copy.number"), with = FALSE]
  template_gr <- GenomicRanges::makeGRangesFromDataFrame(template_data, keep.extra.columns = FALSE)

  NewFormat <- list()
  NewFormat <- lapply(seq_along(Bin_CN_list), function(i) {
    # Extract current Bin_CN_list element
    cell_data <- Bin_CN_list[[i]]

    # Copy the template GRanges object and add "copy.number"
    Bins <- template_gr
    S4Vectors::mcols(Bins)$copy.number <- cell_data$copy.number

    # Compute breakpoints
    shifted_cn <- c(NA, cell_data$copy.number[-nrow(cell_data)])  # 向前平移
    breakpoint_rows <- cell_data$copy.number != shifted_cn & !is.na(shifted_cn)
    breakpoints <- cell_data[breakpoint_rows, c("seqnames", "start", "end", "copy.number"), with = FALSE]


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

  on.exit(future::plan(future::sequential), add = TRUE)

  return(NewFormat)
}


# 1.1: get the copy number matrix
#' Generate copy number matrix for selected cells
#'
#' This function extracts copy number variations from a list of `GRanges` objects
#' and organizes them into an integer matrix. The matrix contains selected cells as columns,
#' with genomic regions (fixed bins) as rows.
#'
#' @param input A named list where each element represents a single cell as a `GRanges` object.
#' @param Template A character vector containing the `cellID`s of selected cells to be included in the matrix.
#'
#' @return An integer matrix:
#'   - Columns represent the selected `cellID`s.
#'   - Rows represent genomic regions, separated into fixed bins.
#' @export
#'
#' @examples
#'
#' file_path <- system.file("extdata", "example_data.rds", package = "cnvTree")
#' Example_data <- changeFormat(file = file_path, core = 4)
#' CNmatrix <- CN_seq(input = Example_data, Template = names(Example_data)[1:10])
#'
#'
CN_seq <- function(input, Template){
    c <- list()
    c <- lapply(Template, function(i) {
      # make sure input[[i]] exist, and include "bins" and "copy.numer" information
      if (!is.null(input[[i]]) && "bins" %in% names(input[[i]]) && "copy.number" %in% names(S4Vectors::mcols(input[[i]][["bins"]]))) {
        as.integer(S4Vectors::mcols(input[[i]][["bins"]])$copy.number)
      } else {
        warning(sprintf("Missing data for template '%s'", i))
        NA_integer_  # 如果數據缺失，用 NA 填充
      }
    })

    c <- as.data.frame(c, stringsAsFactors = FALSE)  # 避免字符串變為因子
    colnames(c) <- Template

    return(c)
  }


