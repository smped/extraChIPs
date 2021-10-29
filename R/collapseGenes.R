#' @title Collapse a vector of gene names
#'
#' @description Collapse a vector of gene names
#'
#' @details
#' Convenience function to collapse a vector of gene names into a character/glue
#' object of length 1. By default, symbols are deduplicated, sorted
#' alpha-numerically and italicised with an underscore.
#'
#' @param x character vector representing gene names
#' @param sort logical(1) Should the names be sorted alphabetically
#' @param dedup logical(1) Should duplicate names be removed
#' @param format character string for markdown formatting as italics/bold
#' @param sep separator between vector elements
#' @param last character string to plast before the last element
#' @param numeric logical(1) sort digits numerically, instead of as strings
#' @param width The maximum width of the string before truncating to ...
#' @param ... passed to \link[stringr]{str_sort}
#'
#' @return a glue object
#'
#' @examples
#' genes <- c("FOXP3", "BRCA1", "TP53")
#' collapseGenes(genes)
#'
#' @importFrom glue glue_collapse
#' @importFrom stringr str_sort
#' @export
collapseGenes <- function(
  x, sort = TRUE, dedup = TRUE, format = "_", sep = ", ", last = " and ",
  numeric = TRUE, width = Inf, ...
) {
  if (length(x) == 0) return(glue(x))
  if (dedup) x <- unique(x)
  if (sort) x <- str_sort(x, numeric = numeric, ...)
  x <- paste0(format, x, format)
  glue_collapse(x, sep = sep, last = last, width = width)
}
