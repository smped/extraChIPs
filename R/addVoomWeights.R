#' @title Add Voom Weights to a SummarizedExperiment
#'
#' @description Add Voom Weights to a Summarized Experiment from counts or CPM
#'
#' @details
#' This adds an assay titled "voomweights" to the SummarizedExperiment object.
#' These can be calculated directly from the "counts" assay or from any other.
#' In particular, this is intended to be used on CPM or normalised CPM values.
#' Calling this directly on the counts assay will use the standard
#' \link[limma]{voom} approach, whilst other assays will assume CPM values
#'
#' @param x A SummarizedExperiment object
#' @param assay The assay to use for estimation of Voom Precision Weights
#' @param design The design matrix
#' @param w0 Vector of initial weights as calculated by
#' \link[limma]{arrayWeights}
#' @param isLogCPM logical(1). Ignored if assay = "counts"
#' @param ... Passed to \link[limma]{voom} and \link[limma]{arrayWeights}
#' @param span See \link[limma]{voom}
#'
#' @return
#' A SummarizedExperiment object
#'
#'
#' @name addVoomWeights
#' @rdname addVoomWeights-methods
#' @export
#'
setGeneric(
  "addVoomWeights",
  function(x, ...){standardGeneric("addVoomWeights")}
)
#' @importFrom SummarizedExperiment assay 'assay<-'
#' @importFrom SummarizedExperiment 'colData<-' 'rowData<-'
#' @importFrom limma voom arrayWeights
#' @rdname addVoomWeights-methods
#' @export
setMethod(
  "addVoomWeights",
  signature = signature(x = "SummarizedExperiment"),
  function(
    x, assay = "counts", design = NULL, w0, isLogCPM = TRUE, ..., span = 0.5
  ) {

    ## NULL specs for design will be handled by downstream functions
    if (!is.null(design)) stopifnot(nrow(design) == ncol(x))
    if (!missing(w0)) stopifnot(length(w0) == ncol(x))
    if (is.null(x$totals)) stop(
      "No library sizes are included in 'x'.",
      "Please add these as the column 'totals'"
    )

    ## Separate the arguments for voom & arrayWeights
    dots <- list(...)
    v_allowed <- setdiff(names(formals(voom)), "save.plot")
    aw_allowed <- names(formals(arrayWeights))

    mat <- assay(x, assay)
    if (assay == "counts" & missing(w0)) {
      voom_args <- list(
        counts = mat, design = design, lib.size = x$totals, span = span,
        save.plot = TRUE
      )
      voom_args <- c(voom_args, dots[intersect(names(dots), v_allowed)])
      v <- do.call(voom, voom_args)
    }
    if (assay == "counts" & !missing(w0)) {

      voom_args <- list(
        counts = mat, design = design, lib.size = x$totals, span = span,
        weights = w0, save.plot = TRUE
      )
      voom_args <- c(voom_args, dots[intersect(names(dots), v_allowed)])
      v <- do.call(voom, voom_args)

      aw_args <- list(object = v, design = design)
      aw_args <- c(aw_args, dots[intersect(names(dots), aw_allowed)])
      aw <- do.call(arrayWeights, aw_args)
      v$weights <- t(aw * t(v$weights))
      v$targets$sample.weights <- aw
    }
    if (assay != "counts"){
      v <- voomWeightsFromCPM(
        mat, design = design, w0 = w0, lib.size = x$totals, isLogCPM = isLogCPM,
        span = span, ...
      )
    }

    assay(x, "voomweights") <- v$weights
    if (!missing(w0)) colData(x)$sample.weights <- v$targets$sample.weights
    rowData(x)$sqrt.stdev <- v$voom.xy$y
    x

  }
)
