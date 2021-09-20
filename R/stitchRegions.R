#' @title Stitch GRanges Within A Given Distance
#'
#' @description Stitch together ranges within a given distance
#'
#' @details
#' Stitches together ranges within a given distance, using any ranges provided
#' for exclusion as barriers between stitched ranges. This may be particularly
#' useful if wanting to stitch enhancers whilst excluding promoters.
#'
#' All inputs and outputs are Genomic Ranges objects
#'
#' @param x Ranges to be stitched together
#' @param exclude Ranges to exclude
#' @param maxgap The maximum distance between ranges to be stitched
#' @param ignore.strand logical
#'
#' @return
#' A GRanges object
#'
#' @examples
#' x <- GRanges(c("chr1:1-10", "chr1:101-110", "chr1:1001-1010", "chr2:1-10"))
#' stitchRanges(x, exclude = GRanges("chr1:200:+"))
#' stitchRanges(x)
#'
#' @importFrom plyranges join_nearest_upstream
#' @importFrom forcats fct_explicit_na
#' @importFrom S4Vectors splitAsList
#' @importFrom GenomicRanges `strand<-`
#' @export
stitchRanges <- function(x, exclude, maxgap = 12500L, ignore.strand = TRUE) {

  ## Argument checks
  stopifnot(is(x, "GRanges"))
  if (missing(exclude)) exclude <- GRanges()
  stopifnot(is(exclude, "GRanges"))
  stopifnot(is.numeric(maxgap))
  stopifnot(is.logical(ignore.strand))

  if (ignore.strand) {
    strand(x) <- "*"
    strand(exclude) <- "*"
  }

  exclude$id <- as.factor(seq_along(exclude))
  out <- join_nearest_upstream(x, exclude)
  ## Join nearest will drop any ranges beyond which there is no exclusionary
  ## range. Add these back & set as End-Of-Chromosome so they're not dropped
  out <- c(out, GenomicRanges::setdiff(x, out))
  out$id <- fct_explicit_na(out$id, "EOC")
  out <- splitAsList(out, f = out$id)
  out <- GenomicRanges::reduce(out, min.gapwidth = maxgap)
  out <- unlist(out)
  names(out) <- c()

  out

}
