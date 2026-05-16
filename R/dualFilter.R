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
#' @param restrict Restrict the calculations to a subset of chromosomes. Will
#' only apply when preparing counts for passing to the internal call to
#' \link[csaw]{scaleControlFilter}
#' @param BPPARAM Settings for running in parallel
#' @param verbose Receive progress messages as different steps of the function
#' are begun
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
#' @importFrom IRanges overlapsAny
#' @importFrom methods is
#' @importFrom csaw windowCounts readParam scaleControlFilter
#' @importFrom csaw filterWindowsProportion reform getWidths
#' @importFrom BiocParallel bpparam bplapply bpisup bpstart bpstop
#' @importFrom edgeR cpm
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom stats quantile
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @importMethodsFrom SummarizedExperiment rowRanges rowData 'rowData<-'
#' @importMethodsFrom SummarizedExperiment colData 'colData<-'
#' @importMethodsFrom SummarizedExperiment assay 'assay<-'
#'
#' @export
dualFilter <- function(
        x, bg = NULL, ref, q = 0.75, logCPM = TRUE, keep.totals = TRUE,
        bin.size = 1e4, prior.count = 2, BPPARAM = bpparam(), restrict = NULL,
        verbose = FALSE
) {

    if (!(keep.totals)) {
        msp <- paste(
            "The 'keep.totals' argument in 'dualFilter()' is deprecated and will",
            "be removed in the next release cycle",
            "Library sizes will then be taken from the complete BAM-level counts,",
            "not the retained windows, which is the more appropriate",
            "behaviour for ChIP-seq library size normalisation.",
            "If you require custom library sizes, set 'x$totals' manually after",
            "calling 'dualFilter()'."
        )
        .Deprecated(msg = msg)
    }

    ## Argument checks
    stopifnot(is(x, "RangedSummarizedExperiment"))
    stopifnot(is(ref, "GRanges"))
    stopifnot(q <= 1, q > 0)
    stopifnot(is.logical(logCPM))

    ## Check the BamFiles exist
    stopifnot("bam.files" %in% colnames(colData(x)))
    stopifnot(all(file.exists(colData(x)$bam.files)))
    bfl <- BamFileList(colData(x)$bam.files)
    names(bfl) <- colnames(x)

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
        stopifnot(all(file.exists(colData(bg)$bam.files)))
        bg_bfl <- BamFileList(colData(bg)$bam.files)
        names(bg_bfl) <- colnames(bg)

        if (is.null(bin.size)) {
            bin.size <- max(20 * max(width(x)), 2e3)
            message("Setting bin.size as ", bin.size)
        }

        rp <- readParam()
        if (!is.null(metadata(x)$param)) rp <- metadata(x)$param
        if (!is.null(restrict)) {
            lv <- slot(rp, "restrict")
            restrict <- intersect(restrict, lv)
            if (length(restrict)) rp <- reform(rp, restrict = restrict)
        }

        if (verbose) message("Reading signal counts...")
        signal_counts <- windowCounts(
            bam.files = bfl,
            spacing = bin.size, filter = 0, param = rp, BPPARAM = BPPARAM
        )
        if (keep.totals) signal_counts$totals <- x[,names(bfl)]$totals
        if (verbose) message("Reading bg counts...")
        bg_counts <- windowCounts(
            bam.files = bg_bfl,
            spacing = bin.size, filter = 0, param = rp, BPPARAM = BPPARAM
        )
        if (keep.totals) bg_counts$totals <- bg[,names(bg_bfl)]$totals
        if (verbose) message("Running scaleControlFilter")
        scf <- scaleControlFilter(signal_counts, bg_counts)

        ## Replicate filterWindowsControl without the copy-on-modify RAM blowout
        bg_totals_scaled <- bg$totals * scf$scale
        ## .scaledAverage only needs the counts assay + totals, not the full RSE
        ## Pass lightweight substitutes and only return the filter.stat
        relative.width <- getWidths(bg) / getWidths(x)
        lib.adjust <- prior.count * mean(bg_totals_scaled) / mean(x$totals)
        if (verbose) message("Running scaledAverage on signal")
        abundances <- .scaledAverage(x, scale = 1, prior.count = prior.count)
        if (verbose) message("Running scaledAverage on background")
        bg.ab <- .scaledAverage(
            bg, scale = relative.width, prior.count = lib.adjust,
            lib.size = bg_totals_scaled
        )  # avoids the totals mutation
        control_filter <- abundances - bg.ab

        ## Now proceed as previously
        q <- sqrt(q)
        cuts$control <- quantile(control_filter[ol], probs = 1 - q)
        keep_control <- control_filter > cuts$control
    }

    ## Apply the filter using the expression percentile. This is quick already
    if (verbose) message("Running filterWindowsProportion")
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


#' @keywords internal
.scaledAverage <- function(y, scale = 1, prior.count = NULL, lib.size = NULL) {
    ## This directly lifts the scaledAverage function from csaw, but allows
    ## for passing of the library size which can help avoid copy-on-modify
    ## issues in functions which call this as a lower-level function
    ## Credit to Aaron Lun for the original
    stopifnot("counts" %in% assayNames(y))
    counts <- assay(y, i = "counts", withDimnames = FALSE)
    if (is.null(lib.size)) lib.size <- y$totals
    if (!is.null(y$norm.factors)) {
        lib.size <- lib.size * y$norm.factors
    }
    dispersion <- 0.05
    if (is.null(prior.count)) {
        prior.count <- formals(edgeR::aveLogCPM.DGEList)$prior.count
    }
    is.zero <- scale == 0
    is.neg <- scale < 0
    failed <- is.zero | is.neg
    if (any(failed)) {
        scale[failed] <- 1
    }
    empty <- matrix(0, nrow(counts), ncol(counts))
    ap <- edgeR::addPriorCount(empty, lib.size = lib.size, prior.count = prior.count)
    ap$y <- ap$y * scale
    ap$y <- ap$y + counts
    ave <- edgeR::mglmOneGroup(y = ap$y, offset = ap$offset, dispersion = dispersion)
    ave <- (ave - log(scale) + log(1e+06))/log(2)
    if (any(failed)) {
        if (length(failed) == 1) {
            ave[] <- ifelse(is.zero, -Inf, NA_real_)
        }
        else {
            ave[is.neg] <- NA_real_
            ave[is.zero] <- -Inf
        }
    }
    return(ave)
}
