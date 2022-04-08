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
#' @param simplify logical(1). Attempt to simplify returned columns where
#' possible
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
reduceMC <- function(x, ignore.strand = FALSE, simplify = TRUE, ...) {
    if (!is(x, "GenomicRanges"))
        stop("'x' must be a GenomicRanges object")
    gr <- GenomicRanges::reduce(x, ignore.strand = ignore.strand, ...)
    .mapMcols2Ranges(gr, x, ignore.strand, simplify)
}
