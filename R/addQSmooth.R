#' @title Smooth Quantile Normalise a SummarizedExperiment
#'
#' @description Smooth Quantile Normalise a SummarizedExperiment
#'
#' @details
#' This function adds the assay 'qsmooth' to a SummarizedExperiment or any
#' derivative class, choosing which existing assay to normalise.
#'
#' The qsmooth weights are returned in the metadata element of the final object.
#'
#' @param x a RangedSummarizedExperiment object
#' @param assay The assay to apply SmoothQuantileNormalisation to
#' @param factor A factor used to group the data. See \link[qsmooth]{qsmooth}
#' @param ... Passed to \link[qsmooth]{qsmooth}
#'
#' @return
#' A SummarizedExperiment with additional metadata and the qsmooth assay
#'
#' @importFrom qsmooth qsmooth qsmoothData qsmoothWeights
#' @importFrom S4Vectors 'metadata<-'
#' @importFrom SummarizedExperiment assay 'assay<-'
#' @importClassesFrom qsmooth qsmooth
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @name addQSmooth
#' @rdname addQSmooth-methods
#' @export
#'
setGeneric(
  "addQSmooth",
  function(x, ...){standardGeneric("addQSmooth")}
)
#' @rdname addQSmooth-methods
#' @export
setMethod(
  "addQSmooth",
  signature = signature(x = "SummarizedExperiment"),
  function(x, assay = "counts", factor, ...) {
    if (missing(factor)) stop("The grouping factor must be supplied")
    mat <- assay(x, assay)
    qs <- qsmooth(mat, group_factor = factor, ...)
    assay(x, "qsmooth") <- qsmoothData(qs)
    metadata(x)$qsmoothWeights <- qsmoothWeights(qs)
    x
  }
)

