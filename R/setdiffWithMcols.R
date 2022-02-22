#' @title Perform setdiff operations retaining mcols
#'
#' @description
#' Perform setdiff operations retaining all mcols from the query range
#'
#' @details
#' This extends the methods provided by \link[GenomicRanges]{setdiff}
#'
#' @param x,y GenomicRanges objects
#' @param ignore.strand If set to TRUE, then the strand of x and y is set to
#' "*" prior to any computation.
#' @param ... Not used
#'
#' @return
#' A GRanges object with all mcols returned form the original object.
#' If a range obtained by setdiff maps back to two or more ranges in the
#' original set of Ranges, mcols will be returned as
#' \link[IRanges]{CompressedList} or
#' \link[S4Vectors]{SimpleList} columns
#'
#' @examples
#' x <- GRanges("chr1:1-100:+")
#' x$id <- "range1"
#' y <- GRanges(c("chr1:51-60:+", "chr1:21-30:-"))
#' setdiffWithMcols(x, y)
#' setdiffWithMcols(x, y, ignore.strand = TRUE)
#'
#' @importClassesFrom IRanges CompressedList
#' @importFrom GenomicRanges setdiff findOverlaps
#' @importFrom S4Vectors splitAsList mcols queryHits subjectHits endoapply
#' @importFrom S4Vectors DataFrame SimpleList
#' @importFrom BiocParallel bplapply bpparam
#'
#' @export
setdiffWithMcols <- function(
  x, y, ignore.strand = FALSE,  ...
) {
  ###################################################################
  ## This should be setup as an S4 method for GRanges, GRangesList ##
  ## Maybe GInteractions objects also?                             ##
  ###################################################################
  if (ncol(mcols(x)) == 0) return(setdiff(x, y, ignore.strand))

  gr <- setdiff(x, y, ignore.strand)
  hits <- findOverlaps(gr, x, ignore.strand = ignore.strand)
  i <- queryHits(hits)
  j <- subjectHits(hits)

  ## If mcols only has one column, a vector will be returned so this handles
  ## the multi-column and single-column case
  DF <- DataFrame(mcols(x)[j,])
  DF <- setNames(DF, names(mcols(x)))

  ## Return columns as Compressed/SimpleList objects if the ranges map back
  ## to multiple ranges in the original object
  needs_list <- any(duplicated(i))
  if (needs_list) {
    ## This will work for mcols which are not lists. Need to figure an approach
    ## out which handles columns which are already lists!!!
    DF <- bplapply(
      DF,
      .returnListColumn, i = i, j = j,
      BPPARAM = bpparam()
    )

  }
  mcols(gr) <- DF
  gr

}


.returnListColumn <- function(x, i, j) {

  out <- splitAsList(x[j], f = as.factor(i))
  names(out) <- c()
  if (is(x, "list_OR_List")) {
    ## Restructure so we have a merged list of the same type as the input
    out <- lapply(out, setNames, nm = c())
    out <- lapply(out, unlist)
    out <- as(out, is(x)[[1]])
  }

  ## Now make sure only unique entries are returned, and unlist if possible
  out <- endoapply(out, unique)
  all_unique <- all(vapply(out, function(x) length(x) == 1, logical(1)))
  if (all_unique) out <- unlist(out)
  out

}

