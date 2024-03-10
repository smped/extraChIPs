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
#' @param y_as_both logical(1) If there are any unstranded regions in y, should
#' these be assigned to both strands. If TRUE unstranded regions can be used
#' to partition stranded regions
#' @param ignore.strand If set to TRUE, then the strand of x and y is set to
#' "*" prior to any computation.
#' @param simplify Pass to chopMC and simplify mcols in the output
#' @param suffix Added to any shared column names in the provided objects
#' @param ... Not used
#'
#' @examples
#' x <- GRanges(c("chr1:1-10", "chr1:6-15"))
#' x$id <- paste0("range", seq_along(x))
#' x
#' y <- GRanges(c("chr1:2-5", "chr1:6-12"))
#' y$id <- paste0("range", seq_along(y))
#' y
#' partitionRanges(x, y)
#'
#' @importFrom S4Vectors mcols queryHits subjectHits DataFrame splitAsList
#' @importFrom S4Vectors endoapply mcols<-
#' @import GenomicRanges
#'
#' @export
#' @rdname partitionRanges-methods
setMethod(
    "partitionRanges", c("GRanges", "GRanges"),
    function(
        x, y, y_as_both = TRUE, ignore.strand = FALSE, simplify = TRUE,
        suffix = c(".x", ".y"), ...
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

        if (all(strand(x) == "*")) ignore.strand <- TRUE

        ## If y has overlapping ranges, the partitioning problem doesn't make
        ## sense. Check for overlaps and exit if there are any
        self_ol <- findOverlaps(y, y, ignore.strand = ignore.strand)
        any_self <- any(queryHits(self_ol) != subjectHits(self_ol))
        if (any_self) stop("'y' cannot have any overlapping ranges")

        ## A key issue here is that setting ignore.strand = TRUE will unstrand
        ## the query 'x'. If the subject 'y' is unstranded, but we need to
        ## retain strand information for 'x', the simplest approach is to
        ## duplicate any unstranded regions in 'y' assigning BOTH strands.
        if (!ignore.strand & y_as_both) {
            y_un <- strand(y) == "*"
            y_pos <- y[y_un]
            strand(y_pos) <- "+"
            y_neg <- y[y_un]
            strand(y_neg) <- "-"
            y <- sort(c(y[!y_un], y_pos, y_neg))
        }

        ## First the regions with no overlap. Can't think of a way which doesn't
        ## use lapply. This will be slow for big datasets. Maybe check bplapply?
        grl_x <- splitAsList(x, seq_along(x))
        names(grl_x) <- c()
        sd <- endoapply(grl_x, setdiffMC, y = y, ignore.strand = ignore.strand)
        sd <- unlist(sd)

        ## Now find the regions with overlaps
        ol <- endoapply(grl_x, setdiffMC, sd, ignore.strand = ignore.strand)
        ol <- unlist(ol)
        hits <- findOverlaps(ol, y, ignore.strand = ignore.strand)

        ## Find the intersections using y as the viewpoint
        out <- GRanges(seqinfo = seqinfo(x))
        if (length(hits) > 0) {
            out <- pintersect(
                y[subjectHits(hits)], ol[queryHits(hits)],
                ignore.strand = ignore.strand
            )
            if (ncol(mcols(ol)) > 0) {
                df_ol <- DataFrame(mcols(ol)[queryHits(hits),])
                names(df_ol) <- names(mcols(ol))
                mcols(out) <- cbind(mcols(out), df_ol)
            }
        }
        out <- c(out, sd)
        out <- sort(out, ignore.strand = ignore.strand)
        keep <- names(mcols(out)) %in% c(names(mcols(x)), names(mcols(y)))
        mcols(out) <- mcols(out)[keep]
        if (simplify) out <- chopMC(out, simplify = TRUE)
        out

    }
)
