#' @title Find the best overlap between GRanges
#'
#' @description Find the best overlap between ranges
#'
#' @details
#' This finds the category in the subject GRanges (y) which has the best overlap
#' with the query GRanges (x).
#' The aim is to produce a character vector for best classifying the query
#' GRanges using an external set of features (e.g. promoters, enhancers etc).
#' If the subject (y) is a GRanges object, the values in the specified column
#' will be used as the category.
#' If the subject (y) is a GRangesList, the names of the list will be used to
#' provide the best match
#'
#' @return
#' Character vector the same length as the supplied GRanges object
#'
#' @param x a GRanges object
#' @param y a named GRangesList or GRanges object with mcol as reference
#' category
#' @param var The variable to use as the category. Not required if `y` is a
#' GRangesList
#' @param ignore.strand logical(1) Passed to \link[GenomicRanges]{findOverlaps}
#' @param missing Value to assign to ranges with no overlap
#' @param min_prop Threshold below which overlaps are discarded
#' @param ... Not used
#'
#' @examples
#' gr <- GRanges("chr1:1-10")
#' gr_cat <- GRanges(c("chr1:2-10", "chr1:5-10"))
#' gr_cat$category <- c("a", "b")
#' propOverlap(gr, gr_cat)
#' bestOverlap(gr, gr_cat, var = "category")
#'
#' grl <- splitAsList(gr_cat, gr_cat$category)
#' lapply(grl, function(x) propOverlap(gr, x))
#' bestOverlap(gr, grl)
#'
#' @importFrom S4Vectors mcols splitAsList
#' @rdname bestOverlap-methods
#' @aliases bestOverlap
#' @export
setMethod(
  "bestOverlap",
  signature = signature(x = "GRanges", y = "GRanges"),
  function(
    x, y, var = NULL, ignore.strand = FALSE, missing = NA_character_,
    min_prop = 0.01, ...
  ) {
    cols <- colnames(mcols(y))
    var <- match.arg(var, cols)
    grl <- splitAsList(y, f = mcols(y)[[var]])
    bestOverlap(x, grl, ignore.strand, missing = missing, min_prop = min_prop)
  }
)
#'
#' @importFrom tibble as_tibble tibble
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect all_of
#' @importFrom dplyr arrange distinct left_join filter
#' @rdname bestOverlap-methods
#' @aliases bestOverlap
#' @export
setMethod(
  "bestOverlap",
  signature = signature(x = "GRanges", y = "GRangesList"),
  function(
    x, y, ignore.strand = FALSE, missing = NA_character_, min_prop = 0.01,
    ...
  ) {

    ## x is a GRanges, y is a GRangesList
    if (is.null(names(y))) stop("'y' must be a named GRangesList")

    p <- lapply(y, function(.y) propOverlap(x, .y, ignore.strand))
    tbl <- as_tibble(p)
    tbl[["id"]] <- seq_along(x)
    tbl_long <- pivot_longer(
      tbl, cols = all_of(names(y)), names_to = "y", values_to = "prop"
    )
    id <- prop <- c() # R CMD check error avoidance
    tbl_long <- dplyr::filter(tbl_long, prop >= min_prop)
    tbl_long <- arrange(tbl_long, id, desc(prop))
    tbl_long <- distinct(tbl_long, id, .keep_all = TRUE) # Keep the highest prop
    tbl_long <- arrange(tbl_long, id)

    out <- tibble(id = seq_along(x))
    out <- left_join(out, tbl_long, by = "id")
    out[["y"]][is.na(out[["y"]])] <- missing

    ## Can't see this happening, but put it here anyway
    if (nrow(out) != length(x)) stop("Output length doesn't match input length")

    out[["y"]]

  }
)
#' @rdname bestOverlap-methods
#' @aliases bestOverlap
#' @export
setMethod(
  "bestOverlap", signature = signature(x = "ANY", y = "ANY"),
  function(x, y, ...) .errNotImp(x, y)
)
