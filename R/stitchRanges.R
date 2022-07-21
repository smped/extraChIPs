#' @title Stitch Ranges within a given distance
#'
#' @description
#' Stitch together ranges within a given distance, using excluded ranges as
#' barriers that cannot be crossed
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
#' x <- GRanges(c("chr1:1-10", "chr1:101-110", "chr1:201-210", "chr2:1-10"))
#' y <- GRanges("chr1:200:+")
#' stitchRanges(x, exclude = y, maxgap = 100)
#'
#' @importFrom forcats fct_explicit_na
#' @importFrom S4Vectors splitAsList subjectHits
#' @importFrom GenomicRanges `strand<-` precede resize shift
#' @export
stitchRanges <- function(x, exclude, maxgap = 12500L, ignore.strand = TRUE) {

  ## Argument checks
  stopifnot(is(x, "GRanges"))
  if (missing(exclude)) exclude <- GRanges()
  stopifnot(is(exclude, "GRanges"))
  stopifnot(is.numeric(maxgap))
  stopifnot(is.logical(ignore.strand))

  if (any(overlapsAny(x, exclude))) {
      warning("Ranges provided in 'x' overlap barrier ranges")
      x <- setdiff(x, exclude, ignore.strand = ignore.strand)
  }

  ## Add point ranges to the end of the `exclude` object to ensure no NA values
  ## are returned by precede
  ## Although this should use seqinfo objects, this is more pragmatic and
  ## allows for a missing seqinfo on either object. Perhaps this should be
  ## added as a check at some point. Enforcing strict compatibility is a good
  ## thing
  all_gr <- c(x, exclude)
  all_gr <- range(all_gr, ignore.strand = ignore.strand)
  chr_lim <- resize(
    all_gr, width = 1, fix = "end", ignore.strand = ignore.strand
  )
  chr_lim <- shift(chr_lim, 1)
  exclude <- sort(c(exclude, chr_lim), ignore.strand = ignore.strand)
  hits <- precede(x, exclude, ignore.strand = ignore.strand)

  out <- splitAsList(x, hits)
  out <- GenomicRanges::reduce(out, min.gapwidth = maxgap)
  out <- unlist(out)
  names(out) <- c()
  out

}
