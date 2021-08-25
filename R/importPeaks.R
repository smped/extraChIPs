#' @title Import peaks
#'
#' @description Import peaks in narrowPeak or broadPeak format
#'
#' @details
#' Peaks are imported from either narrowPeak or broadPeak format as
#' GenomicRanges objects.
#'
#' @param x One or more files to be imported
#' @param type The type of peaks to be imported
#' @param seqinfo A seqinfo object to be applied to the GRanges objects
#' @param pruning.mode How to handle conflicts if supplying a seqinfo object.
#' Defaults to `pruning.mode = "coarse"`. Only "coarse" and "error" are
#' implemented. See \link[GenomeInfoDb]{seqinfo}.
#' @param blacklist A set of ranges to be excluded
#' @param sort logical. Should the ranges be sorted during import
#'
#' @return
#' GRanges or GRangesList depending on the length of the supplied files in `x`
#'
#' @importFrom GenomeInfoDb seqinfo seqnames
#' @importFrom IRanges overlapsAny
#' @export
importPeaks <- function(
  x, type = c("narrow", "broad"),
  seqinfo, pruning.mode = c("coarse", "error"),
  blacklist, sort = TRUE
) {

  stopifnot(file.exists(x))
  type <- match.arg(type)
  n <- length(x)
  if (type == "narrow") {
    out <- lapply(
      x,
      .importSingleNarrow,
      seqinfo = seqinfo, blacklist = blacklist, sort = sort
    )
  }

  if (type == "broad") {
    out <- lapply(
      x,
      .importSingleBroad,
      seqinfo = seqinfo, blacklist = blacklist, sort = sort
    )
  }

  out <- GRangesList(out)
  if (n == 1) out <- out[[1]]
  out

}

.importSingleNarrow <- function(x, seqinfo, blacklist, sort) {

  stopifnot(length(x) == 1)
  stopifnot(.isValidNarrow(x))
  stopifnot(is.logical(sort))
  colNames <- c(
    "seqnames", "start", "end", "name", "score", "strand", "signalValue",
    "pValue", "qValue", "peak"
  )
  classes <- c(
    "character", "numeric", "numeric", "character", "numeric", "character",
    rep("numeric", 4)
  )
  df <- read.table(x, sep = "\t", col.names = colNames, colClasses = classes)
  rownames(df) <- df[["name"]]
  df <- df[-4]

  ## Perform the conversion to a GRanges
  if (!missing(seqinfo)) {

    pruning.mode <- match.arg(pruning.mode)
    stopifnot(is(seqinfo, "Seqinfo"))

    ## Fail if there is a mismatch between any supplied seqnames and seqinfo
    if (pruning.mode == "error")
      stopifnot(df[["seqnames"]] %in% seqnames(seqinfo))

    ## Subset to remove any ranges which are not included in seqinfo
    if (pruning.mode == "coarse")
      df <- subset(df, seqnames %in% seqnames(seqinfo))
    if (nrow(df) == 0) message("No ranges match the supplied seqinfo object")

    gr <- makeGRangesFromDataFrame(
      df, keep.extra.columns = TRUE, seqinfo = seqinfo,
      starts.in.df.are.0based = TRUE
    )

  } else {

    gr <- makeGRangesFromDataFrame(
      df, keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE
    )

  }

  ## Apply the blacklist if supplied
  if (!missing(blacklist)) {
    stopifnot(is(blackList, "GRanges"))
    gr <- gr[!overlapsAny(gr, blacklist)]
  }

  ## Sort if required
  if (sort) gr <- sort(gr)

  gr

}

.importSingleBroad <- function(x, seqinfo, blacklist, sort) {

  stopifnot(length(x) == 1)
  stopifnot(.isValidBroad(x))
  stopifnot(is.logical(sort))
  colNames <- c(
    "seqnames", "start", "end", "name", "score", "strand", "signalValue",
    "pValue", "qValue"
  )
  classes <- c(
    "character", "numeric", "numeric", "character", "numeric", "character",
    rep("numeric", 3)
  )
  df <- read.table(x, sep = "\t", col.names = colNames, colClasses = classes)
  rownames(df) <- df[["name"]]
  df <- df[-4]

  ## Perform the conversion to a GRanges
  if (!missing(seqinfo)) {

    pruning.mode <- match.arg(pruning.mode)
    stopifnot(is(seqinfo, "Seqinfo"))

    ## Fail if there is a mismatch between any supplied seqnames and seqinfo
    if (pruning.mode == "error")
      stopifnot(df[["seqnames"]] %in% seqnames(seqinfo))

    ## Subset to remove any ranges which are not included in seqinfo
    if (pruning.mode == "coarse")
      df <- subset(df, seqnames %in% seqnames(seqinfo))
    if (nrow(df) == 0) message("No ranges match the supplied seqinfo object")

    gr <- makeGRangesFromDataFrame(
      df, keep.extra.columns = TRUE, seqinfo = seqinfo,
      starts.in.df.are.0based = TRUE
    )

  } else {

    gr <- makeGRangesFromDataFrame(
      df, keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE
    )

  }

  ## Apply the blacklist if supplied
  if (!missing(blacklist)) {
    stopifnot(is(blackList, "GRanges"))
    gr <- gr[!overlapsAny(gr, blacklist)]
  }

  ## Sort if required
  if (sort) gr <- sort(gr)

  gr

}

.isValidNarrow <- function(x) {
  r1 <- read.table(x, sep = "\t", nrows = 1)
  nCols <- ncol(r1) == 10
  allNumerics <- suppressWarnings(
    !anyNA(
      vapply(r1[c(2, 3, 5, 7, 8, 9, 10)], as.numeric, numeric(1))
    )
  )
  strandOK <- r1[[6]] %in% c("+", "-", ".")
  all(nCols, allNumerics, strandOK)
}

.isValidBroad <- function(x) {
  r1 <- read.table(x, sep = "\t", nrows = 1)
  nCols <- ncol(r1) == 9
  allNumerics <- suppressWarnings(
    !anyNA(
      vapply(r1[c(2, 3, 5, 7, 8, 9)], as.numeric, numeric(1))
    )
  )
  strandOK <- r1[[6]] %in% c("+", "-", ".")
  all(nCols, allNumerics, strandOK)
}
