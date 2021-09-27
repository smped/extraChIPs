#' @title Smooth Quantile Normalise a SummarizedExperiment
#'
#' @description Smooth Quantile Normalise a SummarizedExperiment
#'
#' @details
#' This function adds the assay 'qsmooth' to a SummarizedExperiment or any
#' derivative class, choosing which existing assay to normalise.
#'
#' A subsample up to 10,000 of the qsmooth weights are returned in the metadata
#' element of the final object.
#'
#' @param x a RangedSummarizedExperiment object
#' @param assay The assay to apply SmoothQuantileNormalisation to
#' @param factor The column used to group the data. Must be one of the columns
#' in the colData element of x
#' @param n_w The maximum number of weights to return in the metadata
#' @param ... Passed to \link[qsmooth]{qsmooth}
#'
#' @return
#' A SummarizedExperiment with additional metadata and the qsmooth assay
#'
#' @examples
#' dat <- cbind(
#'   matrix(rnorm(1000), nrow=100, ncol=5),
#'   matrix(rnorm(1000, .1, .7), nrow=100, ncol=5)
#'  )
#' colnames(dat) <- c(paste0("A", 1:5), paste0("B", 1:5))
#' df <- DataFrame(group = rep(c("A", "B"), each = 5))
#' se <- SummarizedExperiment(
#'   assays = SimpleList(counts = dat),
#'   colData = df
#' )
#' se2 <- addQSmooth(se, "counts", "group")
#'
#' @importFrom qsmooth qsmooth qsmoothData qsmoothWeights
#' @importFrom S4Vectors 'metadata<-'
#' @importFrom SummarizedExperiment assay 'assay<-' colData
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
  function(x, assay = "counts", factor, n_w = 1e4, ...) {
    if (missing(factor)) stop("The grouping factor must be supplied")
    stopifnot(factor %in% colnames(colData(x)))
    mat <- assay(x, assay)
    f <- colData(x)[[factor]]
    qs <- qsmooth(mat, group_factor = f, ...)
    assay(x, "qsmooth") <- qsmoothData(qs)
    ## Add a subset of the weights to the metadata for easy plotting
    ## The original qsmoothPlotWeights limits this to 10K points
    w_df <- data.frame(
      quantile = seq(0, 1, length.out = length(qsmoothWeights(qs))),
      weight = qsmoothWeights(qs)
    )
    n <- min(nrow(w_df), n_w)
    if (n == n_w) {
      i <- sample.int(nrow(w_df), n, replace = FALSE)
      w_df <- dplyr::arrange(w_df[i,], quantile)
    }
    metadata(x)$qsmoothWeights <- as(w_df, "DataFrame")
    x
  }
)

