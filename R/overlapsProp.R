#' @title Find the proportions of an overlapping range
#'
#' @description Find the proportion of a query reange which overlaps the subject
#'
#' @details
#' This behaves similarly to \link[IRanges]{overlapsAny} except the proportion
#' of the query range which overlaps one or more subject ranges is returned
#' instead of a logical vector
#'
#' @return
#' Numeric vector the same length as x
#'
#' @param x,y A GenomicRanges object
#' @param ignore.strand If set to TRUE, then the strand of x and y is set to
#' "*" prior to any computation.
#' @param ... Not used
#'
#' @examples
#' x <- GRanges("chr1:1-10")
#' y <- GRanges("chr1:1-5")
#' overlapsProp(x, y)
#' overlapsProp(y, x)
#'
#' @importFrom GenomicRanges intersect findOverlaps width pintersect reduce
#' @importFrom S4Vectors queryHits subjectHits splitAsList
#' @export
#' @rdname overlapsProp
#' @aliases overlapsProp
setMethod(
  "overlapsProp", c("GRanges", "GRanges"),
  function(x, y, ignore.strand = FALSE, ...) {

    hits <- findOverlaps(x, y, ignore.strand = ignore.strand)
    gr <- pintersect(
      x[queryHits(hits)], y[subjectHits(hits)]
    )
    w <- width(reduce(splitAsList(gr, queryHits(hits))))
    w <- vapply(w, sum, integer(1))
    out <- rep(0, length(x))
    out[unique(queryHits(hits))] <- w
    out / width(x)

  }
)
#' @export
#' @rdname overlapsProp
#' @aliases overlapsProp
setMethod("overlapsProp", c("ANY", "ANY"), function(x, y, ...) .errNotImp(x, y))
