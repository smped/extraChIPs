#' @title Convert to a tibble
#'
#' @description Convert multiple Genomic objects to tibbles
#'
#' @details
#' Quick and dirty conversion into a tibble.
#'
#' By default, GenomicRanges will be returned with the range as a character
#' column called `range` and all mcols parsed as the remaining columns.
#' Seqinfo information will be lost during coercion.
#'
#' Given that names for ranges are considered as rownames in the mcols element,
#' these can be simply parsed by setting `rownames = "id"` in the call to
#' `as_tibble()`
#'
#' When coercing a DataFrame, any Compressed/SimpleList columns will be coerced
#' to S3 list columns.
#' Any GRanges columns will be returned as a character column, losing any
#' additional mcols from these secondary ranges
#'
#' Defined as an S3 method for consistency with existing tidy methods
#'
#' @param x A Genomic Ranges or DataFrame object
#' @param rangeAsChar Convert any GRanges element to a character vector
#' @param name Name of column to use for ranges. Ignored if rangeAsChar =
#' `FALSE`
#' @param ... Passed to [tibble::as_tibble()]
#'
#'
#' @return
#' A \link[tibble]{tibble}
#'
#' @examples
#' gr <- GRanges("chr1:1-10")
#' gr$p_value <- runif(1)
#' names(gr) <- "range1"
#' gr
#' as_tibble(gr)
#' as_tibble(gr, rownames = "id")
#' as_tibble(mcols(gr))
#' as_tibble(seqinfo(gr))
#'
#' @importFrom tibble as_tibble
#' @importFrom methods as
#' @importFrom vctrs vec_proxy
#' @importFrom dplyr mutate across
#' @importFrom tidyselect everything
#'
#' @name as_tibble
#' @rdname as_tibble
#' @export
as_tibble.DataFrame <- function(x, rangeAsChar = TRUE, ...) {
    if (rangeAsChar) {
        ## Identify any columns which are GRanges & convert to character
        grCol <- vapply(x, is, logical(1), class2 = "GRanges")
        x[grCol] <- lapply(x[grCol], as.character)
    }
    df <- as.data.frame(x)
    ## This step will handle any Compressed/SimpleList objects by converting
    ## to a generic list. Build a check for this line in case mutate starts
    ## dropping rownames and converting to a tibble by default
    df <- dplyr::mutate(df, across(everything(), vec_proxy))
    as_tibble(df, ...)
}
#' @importFrom tibble as_tibble
#' @importFrom tidyselect all_of everything
#' @importFrom dplyr select
#' @rdname as_tibble
#' @export
as_tibble.GenomicRanges <- function(
    x, rangeAsChar = TRUE, name = "range", ...
) {
    if (!rangeAsChar) return(as_tibble(as.data.frame(x), ...))
    if (name %in% names(mcols(x)))
        stop("A column named ", name, " already exists. Please choose another.")
    tbl <- as_tibble(mcols(x), rangeAsChar = TRUE, ...)
    tbl[[name]] <- as.character(x)
    dplyr::select(tbl, all_of(name), everything())
}
#' @importFrom tibble as_tibble
#' @importFrom methods slot slotNames
#' @rdname as_tibble
#' @export
as_tibble.Seqinfo <- function(x, ...) {
    nm <- slotNames(x)
    df <- lapply(nm, function(i) slot(x, i))
    names(df) <- nm
    as_tibble(df)
}
