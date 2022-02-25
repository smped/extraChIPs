#' @title Reduce ranges retaining mcols
#'
#' @description Reduce ranges retaining mcols
#'
#' @details
#' This function extends \link[GenomicRanges]{reduce} so that **all** `mcols`
#' are returned in the output.
#' Where the reduced ranges map to multiple ranges in the original range,
#' `mcols` will be returned as `CompressedList` columns.
#'
#' If `simplify = TRUE` columns will be returned as vectors where possible.
#'
#' @param x A GenomicRanges object
#' @param ignore.strand If set to TRUE, then the strand of x and y is set to
#' "*" prior to any computation.
#' @param simplify logical(1). Attempt to simplify returned columns where possible
#' @param ... Passed to \link[GenomicRanges]{reduce}
#'
#' @return
#' A GRanges object
#'
#' @examples
#' x <- GRanges(c("chr1:1-10:+", "chr1:6-12:-"))
#' x$id <- c("range1", "range2")
#' reduceMC(x)
#' reduceMC(x, ignore.strand = TRUE)
#'
#' @importFrom GenomicRanges reduce findOverlaps
#' @importFrom S4Vectors mcols queryHits subjectHits
#' @importFrom S4Vectors DataFrame
#' @importFrom BiocParallel bplapply bpparam
#'
#'
#' @export
#' @rdname reduceMC
#' @aliases reduceMC
setMethod(
  "reduceMC", "GRanges",
  function(x, ignore.strand, simplify, ...) {

    gr <- GenomicRanges::reduce(x, ignore.strand = ignore.strand, ...)
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
#' @rdname reduceMC
#' @aliases reduceMC
setMethod("reduceMC", "ANY", function(x, ...) .errNotImp(x))
