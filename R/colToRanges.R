#' @title Coerce a column to a GRanges object
#'
#' @description Coerce a column to a GRanges object
#'
#' @details
#' Take a data.frame-like object and coerce one column to a GRanges object,
#' setting the remainder as the `mcols`.
#' A particularly useful application of this is when you have a GRanges object
#' with one mcol being a secondary GRanges object.
#'
#' Alternatively, if you have a data.frame with GRanges represented as a
#' character column, this provides a simple method of coercion.
#' In this case, no Seqinfo element will be applied to the GRanges element.
#'
#' @return
#' A GenomicRanges object
#'
#' @examples
#' x <- GRanges(c("chr1:1-10", "chr1:6-15", "chr1:51-60"))
#' df <- data.frame(logFC = rnorm(3), logCPM = rnorm(3,8), p = 10^-rexp(3))
#' gr <- mergeByCol(x, df, col = "logCPM", pval = "p")
#' colToRanges(gr, "keyval_range")
#'
#' @param x A data-frame or GRanges object containing the column to coerce
#' @param var The name of the column to coerce
#' @param ... Not used
#'
#' @importFrom methods as
#' @importFrom GenomicRanges 'mcols<-'
#' @rdname colToRanges-methods
#' @aliases colToRanges
#' @export
setMethod(
  "colToRanges",
  signature = signature(x = "DataFrame", var= "character"),
  function(x, var, ...) {
    stopifnot(var %in% colnames(x))
    gr <- as(x[[var]], "GRanges")
    keep <- setdiff(colnames(x), var)
    mcols(gr) <- x[keep]
    gr
  }
)
#' @importFrom GenomicRanges mcols
#' @rdname colToRanges-methods
#' @aliases colToRanges
#' @export
setMethod(
  "colToRanges",
  signature = signature(x = "GRanges", var= "character"),
  function(x, var, ...) {
    df <- mcols(x)
    colToRanges(df, var, ...)
  }
)
#' @importFrom GenomicRanges mcols
#' @rdname colToRanges-methods
#' @aliases colToRanges
#' @export
setMethod(
  "colToRanges",
  signature = signature(x = "data.frame", var= "character"),
  function(x, var, ...) {
    stopifnot(var %in% colnames(x))
    gr <- as(x[[var]], "GRanges")
    keep <- setdiff(colnames(x), var)
    mcols(gr) <- x[keep]
    gr
  }
)
#' @rdname colToRanges-methods
#' @aliases colToRanges
#' @export
setMethod(
  "colToRanges",
  signature = signature(x = "ANY", var= "ANY"),
  function(x, var, ...) .errNotImp(x, var)
)
