#' @title Convert GRanges to a tibble
#'
#' @description Convert GRanges objects to tibble objects
#'
#' @details
#' Simple conversion of a GRanges object to a tibble. The ranges will be
#' returned as a character column called `ranges`. Seqinfo information will be
#' lost.
#'
#' @param x A Genomic Ranges object
#' @param var Name of column to use for ranges
#' @param drop.existing logical(1). Overwrite any existing columns with the
#' same name as supplied in var
#'
#' @importFrom tibble as_tibble
#' @importFrom GenomicRanges granges
#' @importFrom S4Vectors mcols
#'
#' @return
#' A \link[tibble]{tibble}
#'
#' @examples
#' gr <- GRanges("chr1:1-10")
#' gr$p_value <- runif(1)
#' rangesAsTibble(gr)
#' mcolsAsTibble(gr)
#'
#' @export
rangesAsTibble <- function(x, var = "range", drop.existing = TRUE) {
    stopifnot(is(x, "GRanges"))
    df <- as.data.frame(mcols(x))
    cols <- colnames(df)
    if (!drop.existing & var  %in% cols)
        stop("The column ", var, " already exists")
    df[[var]] <- as.character(granges(x))
    as_tibble(df[unique(c(var, cols))])
}

#' @importFrom tibble as_tibble
#' @importFrom GenomicRanges granges
#' @importFrom S4Vectors mcols
#' @importFrom methods slotNames
#'
#' @export
#' @rdname rangesAsTibble
mcolsAsTibble <- function(x) {
    stopifnot(is(x, "Vector"))
    stopifnot("elementMetadata" %in% slotNames(x))
    df <- as.data.frame(mcols(x))
    as_tibble(df)
}
