#' @title Perform set operations retaining mcols
#'
#' @description
#' Perform set operations retaining all mcols from the query range
#'
#' @details
#' This extends the methods provided by \link[GenomicRanges]{setdiff}
#'
#' @param x,y GenomicRanges objects
#' @param ignore.strand If set to TRUE, then the strand of x and y is set to
#' "*" prior to any computation.
#' @param simplify logical(1) If TRUE, any List columns will be returned as
#' vectors where possible. This can only occur if single, unique entries are
#' present in all initial elements.
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
#' setdiffMC(x, y)
#' setdiffMC(x, y, ignore.strand = TRUE)
#'
#' @importClassesFrom IRanges CompressedList
#' @importFrom GenomicRanges setdiff findOverlaps
#' @importFrom S4Vectors splitAsList mcols queryHits subjectHits endoapply
#' @importFrom S4Vectors DataFrame SimpleList
#' @importFrom BiocParallel bplapply bpparam
#'
#' @export
#' @rdname setoptsMC
#' @aliases setdiffMC
setMethod(
  "setdiffMC", c("GRanges", "GRanges"),
  function(x, y, ignore.strand = FALSE, simplify = TRUE, ...) {
  ###################################################################
    ## This should be setup as an S4 method for GRanges, GRangesList ##
    ## Maybe GInteractions objects also?                             ##
    ###################################################################
    ###################################################################
    ## Write for all setopts: union, intersect, setdiff & reduce
    ###################################################################
    gr <- GenomicRanges::setdiff(x, y, ignore.strand)
    if (ncol(mcols(x)) == 0) return(gr)

    hits <- findOverlaps(gr, x, ignore.strand = ignore.strand)
    i <- queryHits(hits)
    j <- subjectHits(hits)

    ## If mcols only has one column, a vector will be returned so this handles
    ## the multi-column and single-column case
    DF <- DataFrame(mcols(x)[j,])
    DF <- setNames(DF, names(mcols(x)))

    ## Return columns as CompressedList objects if the new ranges map
    ## to multiple ranges in the original object
    if (any(duplicated(i))) {
      DF <- bplapply(
        DF,
        .returnListColumn, i = i, j = j, .simplify = simplify,
        BPPARAM = bpparam()
      )

    }
    mcols(gr) <- DF
    gr

  }
)
#' @export
#' @rdname setoptsMC
#' @aliases setdiffMC
setMethod(
  "setdiffMC", c("ANY", "ANY"),
  function(x, y, ...) .errNotImp(x, y)
)


.returnListColumn <- function(x, i, j, .simplify) {

  out <- splitAsList(x[j], f = as.factor(i))
  names(out) <- c()
  if (is(x, "list_OR_List")) {
    ## Restructure so we have a merged list of the same type as the input
    out <- lapply(out, setNames, nm = c())
    out <- lapply(out, unlist)
    out <- as(out, is(x)[[1]])
  }

  if (.simplify) {
    ## Now make sure only unique entries are returned, and unlist if possible
    out <- endoapply(out, unique)
    all_unique <- all(vapply(out, function(x) length(x) == 1, logical(1)))
    if (all_unique) out <- unlist(out)
  }

  out

}

