#' @title Convert GRanges or DFrame to a tibble
#'
#' @description Convert GRanges or DataFrame objects to tibble objects
#'
#' @details
#' Quick and dirty conversion into a tibble. The ranges will be
#' returned as a character column called `range`. Seqinfo information will be
#' lost.
#'
#' @param x A Genomic Ranges or DataFrame object
#' @param name Name of column to use for ranges. Ignored if rangeAsChar =
#' `FALSE`
#' @param rangeAsChar Convert the GRanges element to a character vector
#' @param ... Not used
#'
#'
#' @return
#' A \link[tibble]{tibble}
#'
#' @examples
#' gr <- GRanges("chr1:1-10")
#' gr$p_value <- runif(1)
#' as_tibble(mcols(gr))
#' as_tibble(gr)
#'
#' @importFrom tibble as_tibble
#' @importFrom methods as
#'
#' @rdname as_tibble-methods
#' @export
as_tibble.DataFrame <- function(x, ...) {
    as_tibble(as.data.frame(x))
}
#' @importFrom rlang ':='
#' @importFrom tibble tibble as_tibble
#' @importFrom dplyr bind_cols
#' @rdname as_tibble-methods
#' @export
as_tibble.GenomicRanges <- function(
    x, name = "range", rangeAsChar = TRUE, ...
) {
    if (!rangeAsChar) return(as_tibble(as.data.frame(x)))
    gr_tbl <- tibble("{name}" := as.character(x))
    if (ncol(mcols(x)) == 0) return(gr_tbl)
    mc_tbl <- as_tibble(mcols(x))
    bind_cols(gr_tbl, mc_tbl)
}
