#' @title Make a set of consensus peaks
#'
#' @description Make a set of consensus peaks based on umber of replicates
#'
#' @details
#' This takes a list of GRanges objects and forms a set of consensus peaks
#' using the minimum proportion of replicates specified
#'
#' @param x A GRangesList
#' @param p The minimum proportion of samples (i.e. elements of `x`) required
#' for a peak to be retained in the output. By default all merged peaks will
#' be returned
#' @param var Additional columns in the `mcols` element to retain
#' @param min.gapwidth,ignore.strand Passed to \link[GenomicRanges]{reduce}
#'
#' @return
#' `GRanges` object with mcols containing a logical vector for every element of
#' x, along with the column `n` which adds all logical columns
#'
#' @examples
#' a <- GRanges("chr1:11-20")
#' a$score <- 1
#' b <- GRanges(c("chr1:18-22", "chr1:1-5"))
#' b$score <- c(0.6, 0.3)
#' grl <- GRangesList(a = a, b = b)
#' makeConsensus(grl)
#' makeConsensus(grl, p = 1)
#' makeConsensus(grl, p = 1, var = "score")
#'
#' @importFrom GenomicRanges granges GRangesList
#' @importFrom IRanges overlapsAny
#' @importFrom S4Vectors "mcols<-" subset mcols endoapply
#' @importFrom methods is
#' @export
makeConsensus <- function(
        x, p = 0, var = NULL, min.gapwidth = 1L, ignore.strand = TRUE
) {

    ## Starting with a GRList
    if (!is(x, "GRangesList")) stop("Input must be a GRangesList")
    if (!is.null(var)) {
        stopifnot(all(var %in% colnames(mcols(x[[1]]))))
        x <- endoapply(
            x,
            function(gr) {
                mcols(gr) <- mcols(gr)[var]
                gr
            }
        )
    } else {
        x <- endoapply(x, granges)
    }

    ## For now, remove all mcols, however reduceMC may be useful if wishing to
    ## retain these for use with plotOverlaps()
    red_ranges <- reduceMC(
        unlist(x), ignore.strand = ignore.strand, min.gapwidth = min.gapwidth
    )
    ol <- lapply(x, function(x) overlapsAny(red_ranges, x))
    ol <- as.data.frame(ol)
    ## Ensure names are strictly retained
    names(ol) <- names(x)
    ol$n <- rowSums(ol)
    mcols(red_ranges) <- cbind(mcols(red_ranges), ol)
    subset(red_ranges, n >= p * length(x))

    ## Now build tests, add to main and change plotOverlaps to use it

}
