#' @title Keep unique ranges and collapse mcols
#'
#' @description Keep unique ranges by 'chopping' mcols
#'
#' @details
#' This function finds unique ranges and chops **all** mcols in a manner similar
#' to \link[tidyr]{chop}.
#' Chopped columns will be returned as `CompressedList` columns, unless
#' `simplify = TRUE` (the default).
#' In this case, columns will be returned as vectors where possible.
#'
#' @param x A GenomicRanges object
#' @param simplify logical(1)
#'
#' @return
#' A GRanges object
#'
#' @examples
#' gr <- GRanges(rep(c("chr1:1-10"), 2))
#' gr$id <- paste0("range", seq_along(gr))
#' gr$gene <- "gene1"
#' gr
#' chopMC(gr)
#'
#'
#' @importFrom tidyr chop all_of
#' @importFrom S4Vectors mcols 'mcols<-'
#' @importFrom GenomeInfoDb seqinfo
#' @importFrom GenomicRanges GRanges
#'
#' @export
chopMC <- function(x, simplify = TRUE) {

    if (!is(x, "GenomicRanges"))
        stop("'x' must be a GenomicRanges object")
    if (ncol(mcols(x)) == 0) return(unique(x))

    tbl <- as_tibble(x, rangeAsChar = TRUE, name = "range")
    tbl <- chop(tbl, -all_of("range"))
    cols <- setdiff(colnames(tbl), "range")
    tbl_list <- as.list(tbl)

    if (simplify) {
        can_simplify <- vapply(
            tbl_list[cols],
            function(x) {
                l <- vapply(x, function(y) length(unique(y)), integer(1))
                all(l == 1)
            },
            logical(1)
        )
        can_simplify <- names(which(can_simplify))
        tbl_list[can_simplify] <- lapply(
            tbl_list[can_simplify],
            function(x) unlist(lapply(x, unique))
        )
    }

    list_cols <- vapply(tbl_list, is, logical(1), class2 = "list")
    tbl_list[list_cols] <- lapply(tbl_list[list_cols], as, "CompressedList")

    gr <- GRanges(tbl_list[["range"]])
    mcols(gr) <- tbl_list[cols]
    gr

  }
