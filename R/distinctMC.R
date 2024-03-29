#' @title Keep distinct ranges and mcols
#'
#' @description Keep distinct ranges by including mcols
#'
#' @details
#'
#' Wrapper to \link[dplyr]{distinct} for `GRanges` objects.
#' Finds unique ranges and mcols in combination and retains only
#' the distinct combinations, in keeping with the `dplyr` function.
#'
#' Will default to `unique(granges(x))` if no columns are provided
#'
#' @param x A GenomicRanges object
#' @param ... `<data-masking>` Passed to \link[dplyr]{distinct}
#' @param .keep_all If `TRUE`, keep all columns in `x`
#'
#' @return
#' A GRanges object
#'
#' @examples
#' gr <- GRanges(rep(c("chr1:1-10"), 2))
#' gr$id <- paste0("range", seq_along(gr))
#' gr$gene <- "gene1"
#' gr
#' distinctMC(gr)
#' distinctMC(gr, gene)
#' distinctMC(gr, gene, .keep_all = TRUE)
#'
#' @importFrom dplyr distinct
#' @importFrom GenomeInfoDb seqinfo
#' @importFrom rlang sym !!
#' @importFrom methods is
#'
#' @export
distinctMC <- function(x, ..., .keep_all = FALSE) {
    stopifnot(is(x, "GRanges"))
    if (length(x) == 0) return(x)
    tbl <- as_tibble(x, rangeAsChar = TRUE)
    tbl <- distinct(tbl, !!sym("range"), ..., .keep_all = .keep_all)
    gr <- colToRanges(tbl, "range", seqinfo = seqinfo(x))
    gr

}
