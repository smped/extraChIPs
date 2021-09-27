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
#' Cutoff values are found such that the `q` of the windows which overlap one of
#' the reference ranges will be returned, along with any others which match the
#' filtering criteria. Cutoff values for both criteria are added to the metadata
#' element of the returned object
#'
#' @param x RangedSummarizedExperiment containing sample counts
#' @param bg RangedSummarizedExperiment containing background/input counts
#' @param ref GRanges object containing ranges where signal is expected
#' @param q The upper percentile of the reference ranges expected to be returned
#' when tuning the filtering criteria
#' @param logCPM logical(1) Add a logCPM assay to the returned data
#' @param BPPARAM Settings for running in parallel
#' @param ... Not used
#'
#' @return
#' A \link[SummarizedExperiment]{RangedSummarizedExperiment} which is a
#' filtered subset of the original object. If requested the assay "logCPM" will
#' be added (`TRUE` by default)
#'
#' @importFrom Rsamtools BamFileList
#' @importFrom BiocIO path
#' @importFrom methods is
#' @importFrom IRanges overlapsAny
#' @importFrom csaw windowCounts readParam scaleControlFilter
#' @importFrom csaw filterWindowsControl filterWindowsProportion
#' @importFrom BiocParallel bpparam
#' @importFrom edgeR cpm
#' @importFrom S4Vectors metadata 'metadata<-'
#' @importFrom stats quantile
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @importMethodsFrom SummarizedExperiment rowRanges colData 'rowData<-'
#' @importMethodsFrom SummarizedExperiment rowData 'rowData<-'
#' @importMethodsFrom SummarizedExperiment assay 'assay<-'
#'
#' @name dualFilter
#' @rdname dualFilter-methods
#' @export
setGeneric("dualFilter", function(x, bg, ref, ...) {
  standardGeneric("dualFilter")
})
#' @rdname dualFilter-methods
#' @export
setMethod(
  "dualFilter",
  signature = signature(
    x = "RangedSummarizedExperiment",
    bg = "RangedSummarizedExperiment",
    ref = "GRanges"
  ),
  function(x, bg, ref, q = 0.99, logCPM = TRUE, BPPARAM = bpparam()) {

    ## Argument checks
    stopifnot(q <= 1, q > 0)
    stopifnot(is.logical(logCPM))
    ## The first two objects must have identical ranges
    stopifnot(all(rowRanges(x) == rowRanges(bg)))

    ## Check the BamFiles exist
    bfl <- BamFileList(colData(x)$bam.files)
    stopifnot(all(file.exists(path(bfl))))
    bg_bfl <- BamFileList(colData(bg)$bam.files)
    stopifnot(all(file.exists(path(bg_bfl))))

    ## Find which ranges overlap the reference
    ol <- overlapsAny(x, ref)
    stopifnot(any(ol))

    rp <- readParam()
    if (!is.null(metadata(x)$param)) rp <- metadata(x)$param

    ## Apply the filter using control samples
    bin_size <- 1e4
    signal_counts <- windowCounts(
      bam.files = bfl,
      spacing = bin_size, filter = 0, param = rp, BPPARAM = BPPARAM
    )
    bg_counts <- windowCounts(
      bam.files = bg_bfl,
      spacing = bin_size, filter = 0, param = rp, BPPARAM = BPPARAM
    )
    scf <- scaleControlFilter(signal_counts, bg_counts)
    control_filter <- filterWindowsControl(
      data = x, background = bg, scale.info = scf
    )$filter
    cuts <- list()
    cuts$control <- quantile(control_filter[ol], probs = 1 - sqrt(q))

    ## Apply the filter using the expression percentile
    prop_filter <- filterWindowsProportion(x)$filter
    cuts$prop <- quantile(prop_filter[ol], probs = 1 - sqrt(q))

    keep <- control_filter > cuts$control & prop_filter > cuts$prop
    out <- x[keep,]
    rowData(out)$overlaps_ref <- ol[keep]
    metadata(out)$cuts <- lapply(cuts, as.numeric)

    if (logCPM){
      assay(out, "logCPM") <- cpm(
        assay(out, "counts"), log = TRUE, lib.size = out$totals
      )
    }

    out

  }
)
