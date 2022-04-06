#' @title Keep distinct ranges and mcols
#'
#' @description Keep distinct ranges by including mcols
#'
#' @details
#' This function finds unique ranges and mcols in combination and retains only
#' the distinct combinations. If no mcols are present defaults to `unique(x)`
#'
#' @param x A GenomicRanges object
#' @param cols The mcols to be used to determining distinct ranges/mcols
#' @param ... Not used
#'
#' @return
#' A GRanges object
#'
#' @examples
#' library(tidyselect)
#' gr <- GRanges(rep(c("chr1:1-10"), 2))
#' gr$id <- paste0("range", seq_along(gr))
#' gr$gene <- "gene1"
#' gr
#' distinctMC(gr)
#' distinctMC(gr, "gene")
#'
#' @importFrom dplyr distinct across
#' @importFrom tidyr all_of
#' @importFrom GenomeInfoDb seqinfo
#'
#' @export
#' @rdname distinctMC-methods
#' @aliases distinctMC
setMethod(
  "distinctMC", "GRanges",
  function(x, cols, ...) {

    if (ncol(mcols(x)) == 0) return(unique(x))
    mc_names <- colnames(mcols(x))
    if (missing(cols)) cols <- mc_names
    sd <- paste(setdiff(cols, mc_names), sep = ", ")
    if (length(sd) > 0)
      stop("Requested columns absent from data: ", sd)

    tbl <- as_tibble(x, name = "range")
    tbl <- distinct(tbl, across(all_of(c("range", cols))))
    gr <- colToRanges(tbl, "range", seqinfo = seqinfo(x))
    gr

  }
)
