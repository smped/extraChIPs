#' @title Partition a set of Genomic Ranges
#'
#' @description Partition a set of Genomic Ranges by another
#'
#' @details
#' The query set of ranges can be broken in regions which stricly overlap
#' a second set of ranges.
#' The complete set of mcols from both initial objects will included in the
#' set of partitioned ranges
#'
#' @return
#' A GRanges object
#'
#' @param x,y GenomicRanges objects
#' @param ignore.strand If set to TRUE, then the strand of x and y is set to
#' "*" prior to any computation.
#' @param suffix Added to any shared column names in the provided objects
#' @param ... Not used
#'
#' @examples
#' x <- GRanges(c("chr1:1-10", "chr1:6-15"))
#' x$id <- paste0("range", seq_along(x))
#' x
#' y <- GRanges(c("chr1:2-5", "chr1:5-12"))
#' y$id <- paste0("range", seq_along(y))
#' y
#' partitionRanges(x, y)
#'
#' @importFrom S4Vectors mcols queryHits subjectHits
#' @importFrom GenomicRanges findOverlaps pintersect setdiff sort
#'
#' @export
#' @rdname partitionRanges-methods
setMethod(
  "partitionRanges", c("GRanges", "GRanges"),
  function(x, y, ignore.strand = FALSE, suffix = c(".x", ".y"), ...) {

    ## First deal with an shared column names in the mcols elements
    cmn_names <- intersect(names(mcols(x)), names(mcols(y)))
    if (length(cmn_names) > 0) {
      i_x <- names(mcols(x)) %in% cmn_names
      names(mcols(x))[i_x] <- paste0(names(mcols(x))[i_x], suffix[[1]])
      i_y <- names(mcols(y)) %in% cmn_names
      names(mcols(y))[i_y] <- paste0(names(mcols(y))[i_y], suffix[[2]])
    }

    ## Now partition the ranges
    hits <- findOverlaps(x, y, ignore.strand = ignore.strand)
    gr <- pintersect(x[queryHits(hits)], y[subjectHits(hits)])
    mcols(gr)[names(mcols(y))] <- mcols(y)[subjectHits(hits),]
    mcols(gr) <- mcols(gr)[colnames(mcols(gr)) != "hit"]

    ## Now add any from x which didn't overlap y
    non_ol <- GenomicRanges::setdiff(x, y, ignore.strand = ignore.strand)
    non_hits <- findOverlaps(non_ol, x)
    non_gr <- pintersect(non_ol[queryHits(non_hits)], x[subjectHits(non_hits)])
    mcols(non_gr) <- mcols(x)[subjectHits(non_hits), ]
    names(mcols(non_gr)) <- names(mcols(x))
    gr <- c(gr, non_gr)
    sort(gr, ignore.strand = ignore.strand)

  }
)
#' @export
#' @rdname partitionRanges-methods
setMethod(
  "partitionRanges", c("ANY", "ANY"), function(x, y, ...) .errNotImp(x, y)
)
