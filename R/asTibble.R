#' @title Convert a GRanges or DataFrame to a tibble
#'
#' @description Convert GRanges or DataFrame objects to tibble objects
#'
#' @details
#' Quick and dirty conversion into a tibble. The ranges will be
#' returned as a character column called `range`. Seqinfo information will be
#' lost.
#'
#' Any Compressed/SimpleList columns will be coerced to S3 list columns
#'
#' Defined as an S3 method for compatibility with existing tidy methods
#'
#' @param x A Genomic Ranges or DataFrame object
#' @param rangeAsChar Convert any GRanges element to a character vector
#' @param name Name of column to use for ranges. Ignored if rangeAsChar =
#' `FALSE`
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
#' @importFrom vctrs vec_proxy
#'
#' @name as_tibble
#' @rdname as_tibble
#' @export
as_tibble.DataFrame <- function(x, rangeAsChar = TRUE, ...) {
    if (rangeAsChar) {
        ## First identify any columns which are GRanges
        ## The convert to character
        grCol <- vapply(x, is, logical(1), class2 = "GRanges")
        x[grCol] <- lapply(x[grCol], as.character)
    }
    df <- as.data.frame(x)
    ## This step will handle any Compressed/SimpleList objects by converting
    ## to a generic list
    df <- lapply(df, vec_proxy)
    as_tibble(df)
}
#' @importFrom rlang ':='
#' @importFrom tibble tibble as_tibble
#' @importFrom dplyr bind_cols
#' @rdname as_tibble
#' @export
as_tibble.GenomicRanges <- function(
    x, rangeAsChar = TRUE, name = "range", ...
) {
    if (name %in% names(mcols(x)))
        stop("A column named ", name, " already exists. Please choose another.")
    if (rangeAsChar) {
        gr <- tibble("{name}" := as.character(x))
    } else {
        gr <- as_tibble(as.data.frame(granges(x)))
    }
    df <- as_tibble(mcols(x), rangeAsChar = rangeAsChar)
    if (nrow(df) == 0) return(gr)
    bind_cols(gr, df)
}
