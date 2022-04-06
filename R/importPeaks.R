#' @title Import peaks
#'
#' @description Import peaks in narrowPeak or broadPeak format
#'
#' @details
#' Peaks are imported from either narrowPeak or broadPeak format as
#' GenomicRanges objects.
#'
#' @param x One or more files to be imported. All files must be of the same
#' type, i.e. narrow or broad
#' @param type The type of peaks to be imported
#' @param blacklist A set of ranges to be excluded
#' @param seqinfo A seqinfo object to be applied to the GRanges objects
#' @param pruning.mode How to handle conflicts if supplying a seqinfo object.
#' Defaults to `pruning.mode = "coarse"`. Only "coarse" and "error" are
#' implemented. See \link[GenomeInfoDb]{seqinfo}.
#' @param sort logical. Should the ranges be sorted during import
#' @param setNames logical Set basename(x) as the name
#' @param ... passed to `sort`
#'
#' @return
#' A GRangesList
#'
#' @examples
#' fl <- system.file(
#' "extdata", "testFiles", "test.narrowPeak", package = "extraChIPs"
#' )
#' gr <- importPeaks(fl, "narrow")
#'
#' @importFrom GenomeInfoDb seqinfo seqnames
#' @importFrom IRanges overlapsAny
#' @importFrom GenomicRanges GRangesList
#' @export
importPeaks <- function(
    x, type = c("narrow", "broad"), blacklist,
    seqinfo, pruning.mode = c("coarse", "error"),
    sort = TRUE, setNames = TRUE, ...
) {

    ## Argument checks
    stopifnot(file.exists(x))
    type <- match.arg(type)
    pruning.mode <- match.arg(pruning.mode)
    stopifnot(is.logical(sort) & is.logical(setNames))
    n <- length(x)
    out <- lapply(
        x,
        .importPeakFile,
        type = type, seqinfo = seqinfo, blacklist = blacklist,
        pruning.mode = pruning.mode
    )

    out <- GRangesList(out)
    if (sort) sort(out, ...)
    if (setNames) names(out) <- basename(x)
    out

}

#' @importFrom GenomicRanges makeGRangesFromDataFrame GRanges
#' @importFrom methods is
#' @importFrom utils read.table
.importPeakFile <- function(x, type, seqinfo, blacklist, pruning.mode) {

    stopifnot(length(x) == 1)
    if (!missing(seqinfo)) {
        stopifnot(is(seqinfo, "Seqinfo"))
    } else {
        seqinfo <- NULL
    }

    ## Handle empty files separately first
    if (file.size(x) == 0) return(GRanges(NULL, seqinfo = seqinfo))

    ## Define the parameters for the two types
    colNames <- c(
        "seqnames", "start", "end", "name", "score", "strand",
        "signalValue", "pValue", "qValue", "peak"
    )
    if (type == "narrow") stopifnot(.isValidNarrow(x))
    if (type == "broad") {
        stopifnot(.isValidBroad(x))
        colNames <- setdiff(colNames, "peak")
    }
    classes <- c("character", "numeric", "numeric", "character", "numeric",
                 "character", rep("numeric", 4))[seq_along(colNames)]

    ## Parse
    df <- read.table(x, sep = "\t", col.names = colNames, colClasses = classes)
    rownames(df) <- df[["name"]]
    df <- df[-4] # The name column has been set as rownames

    ## Deal with the Seqinfo object
    if (!is.null(seqinfo)) {

        ## Fail if there is a mismatch between any seqnames and seqinfo
        if (pruning.mode == "error")
            stopifnot(df[["seqnames"]] %in% seqnames(seqinfo))
        ## Subset to remove any ranges which are not included in seqinfo
        if (pruning.mode == "coarse")
            df <- subset(df, seqnames %in% seqnames(seqinfo))
        if (nrow(df) == 0)
            message("No ranges match the supplied seqinfo object")
    }

    ## Form the object
    gr <- makeGRangesFromDataFrame(
        df, TRUE, seqinfo = seqinfo, starts.in.df.are.0based = TRUE
    )

    ## Apply the blacklist if supplied
    if (missing(blacklist)) blacklist <- GRanges(seqinfo = seqinfo)
    stopifnot(is(blacklist, "GRanges"))
    gr[!overlapsAny(gr, blacklist)]

}

#' @importFrom utils read.table
.isValidNarrow <- function(x) {
    r1 <- read.table(x, sep = "\t", nrows = 1)
    nCols <- ncol(r1) == 10
    if (!nCols) return(FALSE)
    charCols <- c(1, 4)
    allChars <- all(vapply(r1[charCols], is.character, logical(1)))
    numericCols <- c(2, 3, 5, 7:10)
    allNumerics <- all(vapply(r1[numericCols], is.numeric, logical(1)))
    strandOK <- r1[[6]] %in% c("+", "-", ".")
    all(allChars, allNumerics, strandOK)
}

#' @importFrom utils read.table
.isValidBroad <- function(x) {
    r1 <- read.table(x, sep = "\t", nrows = 1)
    nCols <- ncol(r1) == 9
    if (!nCols) return(FALSE)
    charCols <- c(1, 4)
    allChars <- all(vapply(r1[charCols], is.character, logical(1)))
    numericCols <- c(2, 3, 5, 7:9)
    allNumerics <- all(vapply(r1[numericCols], is.numeric, logical(1)))
    strandOK <- r1[[6]] %in% c("+", "-", ".")
    all(allChars, allNumerics, strandOK)
}
