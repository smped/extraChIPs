#' @title Keep distinct ranges and mcols
#'
#' @description Keep distinct ranges by including mcols
#'
#' @details
#' This function finds unique ranges and mcols in combination and retains only
#' the distinct combinations
#'
#' @param x A GenomicRanges object
#' @param .across `<tidy-select>` Passed to \link[dplyr]{across}
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
#' distinctMC(gr, all_of("gene"))
#'
#' @importFrom dplyr distinct across
#' @importFrom tidyr everything all_of
#' @importFrom GenomeInfoDb seqinfo
#'
#' @export
#' @rdname distinctMC-methods
#' @aliases distinctMC
setMethod(
  "distinctMC", "GRanges",
  function(x, .across = everything(), ...) {

    tbl <- as_tibble(x, name = "range")
    tbl <- distinct(tbl, across(c(all_of("range"), .across)))
    gr <- colToRanges(tbl, "range", seqinfo = seqinfo(x))
    gr

  }
)
#' @export
#' @rdname distinctMC-methods
#' @aliases distinctMC
setMethod("distinctMC", "ANY", function(x, ...) .errNotImp(x))
