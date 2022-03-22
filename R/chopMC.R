#' @title Keep unique ranges
#'
#' @description Keep unique ranges by 'chopping' mcols
#'
#' @details
#' This function finds unique ranges and chops all mcols in a manner similar
#' to \link[tidyr]{chop}.
#' Chopped columns will be returned as `CompressedList` columns, unless
#' `simplify = TRUE` (the default).
#' In this case, columns will be returned as vectors where possible.
#'
#' @param x A GenomicRanges object
#' @param simplify logical(1). Attempt to simplify returned columns where possible
#' @param ... Not used
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
#' @importFrom S4Vectors mcols
#' @importFrom GenomeInfoDb seqinfo
#'
#' @export
#' @rdname chopMC-methods
#' @aliases chopMC
setMethod(
  "chopMC", "GRanges",
  function(x, simplify = TRUE, ...) {

    tbl <- as_tibble(x, name = "range")
    tbl <- chop(tbl, -all_of("range"))
    gr <- colToRanges(tbl, "range", seqinfo = seqinfo(x))
    mc <- lapply(
      mcols(gr),
      function(x) {
        un <- lapply(x, function(x) unique(unlist(x)))
        un <- as(un, "CompressedList")
        if (simplify & all(vapply(un, length, integer(1)) == 1)) {
          un <- unlist(un)
        }
        un
      }
    )
    mcols(gr) <- mc
    gr

  }
)
#' @export
#' @rdname chopMC-methods
#' @aliases chopMC
setMethod("chopMC", "ANY", function(x, ...) .errNotImp(x))

