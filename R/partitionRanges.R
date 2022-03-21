#' @title Partition a set of Genomic Ranges
#'
#' @description Partition a set of Genomic Ranges by another
#'
#' @details
#' The query set of ranges can be broken in regions which strictly overlap
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
#' @param simplify logical(1). Simplify any `CompressedList` columns as vectors
#' where possible
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
#' @importFrom S4Vectors mcols queryHits subjectHits DataFrame
#' @importFrom GenomicRanges findOverlaps sort psetdiff pintersect
#'
#' @export
#' @rdname partitionRanges-methods
setMethod(
  "partitionRanges", c("GRanges", "GRanges"),
  function(
    x, y, ignore.strand = FALSE, simplify = TRUE, suffix = c(".x", ".y"), ...
  ) {

    if (length(y) == 0) return(x)

    ## First deal with an shared column names in the mcols elements
    cmn_names <- intersect(names(mcols(x)), names(mcols(y)))
    if (length(cmn_names) > 0) {
      i_x <- names(mcols(x)) %in% cmn_names
      names(mcols(x))[i_x] <- paste0(names(mcols(x))[i_x], suffix[[1]])
      i_y <- names(mcols(y)) %in% cmn_names
      names(mcols(y))[i_y] <- paste0(names(mcols(y))[i_y], suffix[[2]])
    }

    ## If y has overlapping ranges, the partitioning problem doesn't make sense
    ## First run a reduction on y, keeping mcols
    red_y <- reduceMC(y, simplify = TRUE, ignore.strand = ignore.strand)
    ol <- findOverlaps(x, red_y, ignore.strand = ignore.strand)

    ## Find the differences using parallel setdiff
    sd <- psetdiff(
      x[queryHits(ol)], red_y[subjectHits(ol)], ignore.strand = ignore.strand
    )
    mcols(sd) <- DataFrame(mcols(x)[queryHits(ol),])
    colnames(mcols(sd)) <- colnames(mcols(x))
    # sd_ol <- findOverlaps(sd, x, ignore.strand = ignore.strand)
    # sd <- sd[queryHits(sd_ol)]
    # mcols(sd) <- DataFrame(mcols(x)[subjectHits(sd_ol),])

    ## Find the intersection of each range using parallel intersect
    int <- pintersect(
      x[queryHits(ol)], red_y[subjectHits(ol)], ignore.strand = ignore.strand
    )
    mcols(int) <- mcols(int)[colnames(mcols(x))]

    ## Join and map the columns from y
    out <- GenomicRanges::sort(c(sd, int), ignore.strand = ignore.strand)
    df_y <- mcols(.mapMcols2Ranges(out, red_y, ignore.strand, simplify))
    mcols(out) <- cbind(mcols(out), df_y)

    ## Now collapse any identical ranges
    chopRanges(out, simplify = simplify)

  }
)
#' @export
#' @rdname partitionRanges-methods
setMethod(
  "partitionRanges", c("ANY", "ANY"), function(x, y, ...) .errNotImp(x, y)
)

