#' @title Import peaks
#'
#' @description Import peaks in narrowPeak or broadPeak format
#'
#' @details
#' Peaks are imported from either narrowPeak, broadPeak or bed format as
#' individual GenomicRanges objects, before being incorporated as a single
#' GRangesList. If importing bed files, \link[rtracklayer]{import.bed} will be
#' called internally.
#'
#' @param x One or more files to be imported. All files must be of the same
#' type, i.e. narrow or broad
#' @param type The type of peaks to be imported. All files are assumed to be
#' the same format
#' @param blacklist A set of ranges to be excluded
#' @param seqinfo A seqinfo object to be applied to the GRanges objects
#' @param pruning.mode How to handle conflicts if supplying a seqinfo object.
#' Defaults to `pruning.mode = "coarse"`. Only "coarse" and "error" are
#' implemented. See \link[GenomeInfoDb]{seqinfo}.
#' @param sort logical. Should the ranges be sorted during import
#' @param setNames logical Set basename(x) as the name
#' @param centre Add the position of the peak centre. Ignored unless
#' importing narrowPeak files.
#' @param nameRanges Set the values in the 'name' column as the range name.
#' Ignored if parsing a bed file
#' @param ... passed to `sort`
#'
#' @return
#' A GRangesList
#'
#' @examples
#' fl <- system.file(
#'     c("extdata/ER_1.narrowPeak", "extdata/ER_2.narrowPeak"),
#'     package = "extraChIPs"
#' )
#' peaks <- importPeaks(fl)
#' peaks
#'
#' @importFrom GenomeInfoDb seqinfo seqnames
#' @import GenomicRanges
#' @export
importPeaks <- function(
        x, type = c("narrow", "broad", "bed"), blacklist,
        seqinfo, pruning.mode = c("coarse", "error"),
        sort = TRUE, setNames = TRUE, centre = FALSE, nameRanges = TRUE,
        ...
) {

    ## Argument checks
    stopifnot(file.exists(x))
    type <- match.arg(type)
    if (type %in% c("narrow", "broad")) type <- paste0(type, "Peak")
    pruning.mode <- match.arg(pruning.mode)
    logicalArgs <- c(sort, setNames, centre, nameRanges)
    stopifnot(is.logical(logicalArgs))
    out <- lapply(
        x,
        .importPeakFile,
        type = type, seqinfo = seqinfo, blacklist = blacklist,
        pruning.mode = pruning.mode, centre = centre, nameRanges = nameRanges
    )

    out <- GRangesList(out)
    if (sort) sort(out, ...)
    if (setNames) names(out) <- basename(x)
    out

}

#' @import GenomicRanges
#' @importFrom IRanges overlapsAny
#' @importFrom methods is
#' @importFrom utils read.table
#' @importFrom rtracklayer import
.importPeakFile <- function(
        x, type, seqinfo, blacklist, pruning.mode, centre, nameRanges
) {

    stopifnot(length(x) == 1)
    if (!missing(seqinfo)) {
        stopifnot(is(seqinfo, "Seqinfo"))
    } else {
        seqinfo <- NULL
    }

    ## Handle empty files separately first
    if (file.size(x) == 0) return(GRanges(NULL, seqinfo = seqinfo))

    ## Import
    gr <- rtracklayer::import(x, format = type)
    if (!is.null(seqinfo)) seqinfo(gr, pruning.mode = pruning.mode) <- seqinfo
    if (nameRanges & all(!is.na(gr$name))) {
        names(gr) <- gr$name
        gr$name <- NULL
    }
    if (centre & type == "narrowPeak") gr$centre <- start(gr) + gr$peak

    ## Apply the blacklist if supplied
    if (missing(blacklist)) blacklist <- GRanges()
    stopifnot(is(blacklist, "GRanges"))
    gr[!overlapsAny(gr, blacklist)]


}
