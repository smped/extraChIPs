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
#' @importFrom GenomicRanges intersect findOverlaps width
#' @importFrom dplyr left_join group_by summarise
#' @export
#' @rdname overlapsProp
#' @aliases overlapsProp
setMethod(
  "overlapsProp", c("GRanges", "GRanges"),
  function(x, y, ignore.strand = FALSE, ...) {

    gr <- GenomicRanges::intersect(x, y, ignore.strand)
    hits <- as.data.frame(findOverlaps(x, gr, ignore.strand = ignore.strand))
    hits$width <- width(gr)[hits$subjectHits]
    query_df <- data.frame(queryHits = seq_along(x))
    merged_df <- left_join(query_df, hits, by = "queryHits")
    merged_df$width[is.na(merged_df$width)] <- 0
    merged_df <- group_by(merged_df, queryHits)
    merged_df <- summarise(merged_df, width = sum(width))
    merged_df$width / width(x)

  }
)
#' @export
#' @rdname overlapsProp
#' @aliases overlapsProp
setMethod("overlapsProp", c("ANY", "ANY"), function(x, y, ...) .errNotImp(x, y))
