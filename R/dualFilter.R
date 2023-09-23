#' @title Apply two filters to sliding windows
#'
#' @description Apply two filters to counts generated using sliding windows
#'
#' @details
#' This function will take sliding (or tiling) windows for it's input as a
#' \link[SummarizedExperiment]{RangedSummarizedExperiment} object. The dual
#' strategy of applying \link[csaw]{filterWindowsControl} and
#' \link[csaw]{filterWindowsProportion} will then be applied. A set of
#' reference ranges for which signal is expected is used to refine the
#' filtering criteria.
#'
#' Cutoff values are found for both signal relative to input and overall signal,
#' such that the `100*q%` of the (sliding) windows which overlap a reference
#' range will be returned, along with any others which match the
#' dual filtering criteria.
#' In general, higher values of `q` will return more windows as those with weak
#' signal and a marginal overlap with a reference range will be returned.
#' Lower values will ensure that fewer windows, generally with the strongest
#' signal, are retained.
#' Cutoff values for both criteria are added to the metadata
#' element of the returned object.
#'
#' If setting `bg = NULL` the \link[csaw]{filterWindowsControl} step will be
#' ignored and only the \link[csaw]{filterWindowsProportion} will be used.
#' This should only be performed if no Input sample is available.
#'
#' **Please note** that the any `.bam` files referred to in the supplied objects
#' **must** be accessible to this function. It will not run on a separate
#' machine or file structure to that which the original sliding windows were
#' prepared. Please see the example/vignette for runnable code.
#'
#' @param x RangedSummarizedExperiment containing sample counts
#' @param bg RangedSummarizedExperiment containing background/input counts,
#' or alternate method for selecting samples from within x, such as a logical,
#' numeric or character vector
#' @param ref GRanges object containing ranges where signal is expected
#' @param q The upper percentile of the reference ranges expected to be returned
#' when tuning the filtering criteria
#' @param logCPM logical(1) Add a logCPM assay to the returned data
#' @param keep.totals logical(1) Keep the original library sizes or replace
#' using only the retained windows
#' @param bin.size Bin sizes when calling \link[csaw]{filterWindowsControl}.
#' If not specified will default to the largest of 2000bp or 10x the window
#' size
#' @param prior.count Passed to \link[csaw]{filterWindowsControl} and
#' \link[csaw]{filterWindowsProportion}
#' @param BPPARAM Settings for running in parallel
#'
#' @return
#' A \link[SummarizedExperiment]{RangedSummarizedExperiment} which is a
#' filtered subset of the original object. If requested the assay "logCPM" will
#' be added (`TRUE` by default)
#'
#' @examples
#' \donttest{
#' ## Taken from the differential_binding vignette
#' library(tidyverse)
#' library(Rsamtools)
#' library(csaw)
#' library(BiocParallel)
#' library(rtracklayer)
#' ## For this function we need a set of counts using sliding windows and the
#' ## original BamFiles from which they were taken
#' ## First we'll set up the bam file list
#' bfl <- system.file(
#'     "extdata", "bam", c("ex1.bam", "ex2.bam", "input.bam"), package = "extraChIPs"
#'     ) %>%
#'     BamFileList() %>%
#'     setNames(c("ex1", "ex2", "input"))
#'
#' ## Then define the readParam settings for csaw::readParam()
#' rp <- readParam(
#'     pe = "none",
#'     dedup = TRUE,
#'     restrict = "chr10"
#' )
#'
#' ## Now we can form our sliding window object with the counts.
#' wincounts <- windowCounts(
#'     bam.files = bfl,
#'     spacing = 60,
#'     width = 180,
#'     ext = 200,
#'     filter = 1,
#'     param = rp
#' )
#' ## As this is a subset of reads, add the initial library sizes for accuracy
#' ## Note that this step is not normally required
#' wincounts$totals <- c(964076L, 989543L, 1172179L)
#'
#' ## We should also update the metadata for our counts
#' wincounts$sample <- colnames(wincounts)
#' wincounts$treat <- as.factor(c("ctrl", "treat", NA))
#' colData(wincounts)
#'
#' ## The function dualFilter requires a set of peaks which will guide the
#' ## filtering step. This indicate where genuine signal is likely to be found
#' ## and will perform the filtering based on a) signal above the input, and
#' ## b) The overall signal level, using the guide set of peaks to inform the
#' ## cutoff values for inclusion
#' peaks <- import.bed(
#'     system.file("extdata", "peaks.bed.gz", package = "extraChIPs")
#' )
#' filtcounts <- dualFilter(
#'     x = wincounts, bg = "input", ref = peaks,
#'     q = 0.8 # Better to use q = 0.5 on real data
#' )
#' filtcounts
#'
#' }
#'
#' @importFrom Rsamtools BamFileList ScanBamParam countBam
#' @importFrom BiocIO path
#' @importFrom IRanges overlapsAny
#' @importFrom methods is
#' @importFrom csaw windowCounts readParam scaleControlFilter
#' @importFrom csaw filterWindowsControl filterWindowsProportion
#' @importFrom BiocParallel bpparam bplapply bpisup bpstart bpstop
#' @importFrom edgeR cpm
#' @importFrom S4Vectors metadata 'metadata<-'
#' @importFrom stats quantile
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @importMethodsFrom SummarizedExperiment rowRanges rowData 'rowData<-'
#' @importMethodsFrom SummarizedExperiment colData 'colData<-'
#' @importMethodsFrom SummarizedExperiment assay 'assay<-'
#'
#' @export
dualFilter <- function(
    x, bg = NULL, ref, q = 0.5, logCPM = TRUE, keep.totals = TRUE,
    bin.size = NULL, prior.count = 2, BPPARAM = bpparam()
) {

    ## Argument checks
    stopifnot(is(x, "RangedSummarizedExperiment"))
    stopifnot(is(ref, "GRanges"))
    stopifnot(q <= 1, q > 0)
    stopifnot(is.logical(logCPM))

    ## Check the BamFiles exist
    stopifnot("bam.files" %in% colnames(colData(x)))
    bfl <- BamFileList(colData(x)$bam.files)
    names(bfl) <- colnames(x)
    stopifnot(all(file.exists(path(bfl))))

    ## Check BiocParallel is ready to go
    if (!bpisup(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    ## Find which ranges overlap the reference
    ol <- overlapsAny(x, ref)
    stopifnot(any(ol))

    cuts <- list()
    keep_control <- TRUE
    ## Add a routine for when no Input sample is provided
    if (is.null(bg)) {
        msg <- paste(
            "No Input/BG samples provided.",
            "Only the 'filterWindowsProportion' step will be run"
        )
        message(msg)
    } else {

                ## Check the different options for providing the input sample
        if (!is(bg, "RangedSummarizedExperiment")) {
            ## Perform selection by numeric given colnames are optional
            i <- bg # Default if numeric
            if (is.logical(bg)) i <- which(bg)
            if (is.character(bg)) i <- which(colnames(x) %in% bg)
            stopifnot(is.numeric(i) & length(i) > 0)
            bg <- x[,i] # Form bg as a RangedSE subsetting x
            x <- x[,-i]
            bfl <- bfl[-i]
        }
        ## Now the output should always be a RangedSummarizedExperiment
        stopifnot(is(bg, "RangedSummarizedExperiment"))
        stopifnot(all(rowRanges(x) == rowRanges(bg)))
        stopifnot("bam.files" %in% colnames(colData(bg)))
        bg_bfl <- BamFileList(colData(bg)$bam.files)
        stopifnot(all(file.exists(path(bg_bfl))))
        names(bg_bfl) <- colnames(bg)

        ## Apply the filter using control samples using the csaw method which
        ## counts all bam files across the entire genome again
        ## This is currently very slow and an alternative method may be better
        ## Maybe it's worth trying to restrict to a subset of chromosomes as the
        ## steps are 1) Calculate the normalisation factor for bg counts, then
        ## 2) recalculate the filter. Perhaps subsetting the genome will give
        ## similar results. TODO...
        rp <- readParam()
        if (!is.null(metadata(x)$param)) rp <- metadata(x)$param
        if (is.null(bin.size)) bin.size <- max(10 * max(width(x)), 2e3)
        signal_counts <- windowCounts(
            bam.files = bfl,
            spacing = bin.size, filter = 0, param = rp, BPPARAM = BPPARAM
        )
        if (keep.totals) signal_counts$totals <- x[,names(bfl)]$totals
        bg_counts <- windowCounts(
            bam.files = bg_bfl,
            spacing = bin.size, filter = 0, param = rp, BPPARAM = BPPARAM
        )
        if (keep.totals) bg_counts$totals <- bg[,names(bg_bfl)]$totals
        scf <- scaleControlFilter(signal_counts, bg_counts)
        if (!keep.totals) {
            x$totals <- signal_counts$totals
            bg$totals <- bg_counts$totals
        }
        control_filter <- filterWindowsControl(
            data = x, background = bg, prior.count = prior.count,
            scale.info = scf
        )$filter
        q <- sqrt(q)
        cuts$control <- quantile(control_filter[ol], probs = 1 - q)
        keep_control <- control_filter > cuts$control
    }

    ## Apply the filter using the expression percentile. This is quick already
    prop_filter <- filterWindowsProportion(x, prior.count = prior.count)$filter
    cuts$prop <- quantile(prop_filter[ol & keep_control], probs = 1 - q)

    keep <- keep_control & prop_filter > cuts$prop
    out <- x[keep,]
    rowData(out)$overlaps_ref <- ol[keep]
    metadata(out)$cuts <- lapply(cuts, as.numeric)
    colData(out) <- droplevels(colData(out))

    if (!keep.totals) {
        gr <- rowRanges(out)
        totals <- bplapply(bfl, countBam, param = ScanBamParam(which = gr))
        out$totals <- vapply(totals, function(x) sum(x$records), numeric(1))
    }

    if (logCPM) {
        lcpm <- cpm(assay(out, "counts"), log = TRUE, lib.size = out$totals)
        assay(out, "logCPM") <- lcpm
    }

    out

}
