#' @title Estimate Voom Precision Weights Directly From CPM Values
#'
#' @description Estimate Voom Precision Weights Directly From CPM Values
#'
#' @details
#' This function takes CPM or logCPM values and estimates the precision weights
#' as would be done by providing counts directly to the \code{\link{voom}}
#' function.
#' Using this function enables the use of logCPM values which have been
#' normalised using other methods such as Conditional-Quantile or
#' Smooth-Quantile Normalisation.
#'
#' The precision weights are returned as part of the \code{EList} output, and
#' these are automatically passed to the function \code{\link{lmFit}} during
#' model fitting.
#' This will ensure that the mean-variance relationship is appropriate for
#' the linear modelling steps as performed by limma.
#'
#' Initial sample weights can be passed to the function, and should be
#' calculated using \code{\link{arrayWeights}} called on the normalised logCPM
#' values.
#' The returned sample weights will be different to these, given that the
#' function \code{\link{voomWithQualityWeights}} performs two rounds of
#' estimation.
#' The first is on the initial data, with the inappropriate mean-variance
#' relationship, whilst the second round is after incorporation of the precision
#' weights.
#'
#' @param cpm Matrix of CPM or logCPM values
#' @param design The design matrix for the experiment
#' @param w0 Initial vector of sample weights. Should be calculated using
#' \code{\link{arrayWeights}}
#' @param lib.size Initial library sizes. Must be provided as these are no
#' estimable from CPM values
#' @param isLogCPM logical(1). Indicates whether the data is log2 transformed
#' already. Most commonly (e.g. if using the output of cqn) it will be,
#' @param span Width of the smoothing window used for the lowess mean-variance
#' trend. Expressed as a proportion between 0 and 1.
#' @param ... Passed to lmFit internally
#'
#' @return
#' An object of class \code{EList} as would be output by voom.
#' Importantly, there will be no \code{genes} element, although this can be
#' added later.
#' Similarly, the returned \code{targets} element will only contain sample
#' names and library sizes.
#' This can be incorporated with any other metadata as required.
#'
#' Plotting data is always returned, noting the the value \code{sx} has
#' been offset by the library sizes and will be simple logCPM values.
#' As such, the fitted \code{Amean} is also returned in this list element.
#'
#' If initial sample weights were provided, modified weights will also be
#' returned, as the initial function \code{\link{voomWithQualityWeights}}
#' performs two rounds of estimation of sample weights.
#' Here we would simply provide the initial weights a priori, with the
#' second round performed within the function.
#' Importantly, this second round of sample weight estimation uses the precision
#' weights ensuring the correct mean-variance relationship is used for the final
#' estimation of sample weights
#'
#' @examples
#' bamFiles <- system.file("exdata", c("rep1.bam", "rep2.bam"), package="csaw")
#' wc <- csaw::windowCounts(bamFiles, filter=1)
#' cpm <- edgeR::cpm(wc, log = TRUE)
#' el <- voomWeightsFromCPM(cpm, lib.size = wc$totals)
#'
#' @importFrom limma lmFit arrayWeights
#' @importFrom stats approxfun lowess
#' @importFrom methods new
#' @importClassesFrom limma EList
#'
#' @export
voomWeightsFromCPM <- function(
  cpm, design = NULL, w0 = NULL, lib.size = NULL, isLogCPM = TRUE,
  span = 0.5, ...
){

  ## Checks taken from voom internals
  n <- nrow(cpm)
  if (n < 2L)
    stop("Need at least two genes to fit a mean-variance trend")
  m <- min(cpm)
  if (is.na(m))
    stop("NA values not allowed")
  if (m < 0 & !isLogCPM)
    stop("Negative CPM values not allowed")
  if (m == 0 & !isLogCPM)
    stop("Please ensure an offset is used for estimation of CPM values.")

  ## Sort out the design matrix
  if (is.null(design)) {
    design <- matrix(1, ncol(cpm), 1)
    rownames(design) <- colnames(cpm)
    colnames(design) <- "GrandMean"
  }

  ## Library sizes must be supplied & valid
  if (is.null(lib.size))
    stop("Library sizes must be provded as these cannot estimated from CPM")
  i <- ncol(cpm)
  stopifnot(length(lib.size) == i)
  stopifnot(is.numeric(lib.size) & all(lib.size > 0))

  ## Check the initial weights
  if (!is.null(w0))
    stopifnot(length(w0) == i & is.numeric(w0))

  ## Make sure we are on the log scale
  if (!isLogCPM)
    cpm <- log2(cpm)

  ## Now take the main body from voom
  fit <- lmFit(cpm, design, weights = w0, ...)
  if (is.null(fit$Amean))
    fit$Amean <- rowMeans(cpm, na.rm = TRUE)
  sx <- fit$Amean + mean(log2(lib.size + 1)) - log2(1e+06)
  sy <- sqrt(fit$sigma)
  l <- lowess(sx, sy, f = span)
  f <- approxfun(l, rule = 2, ties = list("ordered", mean))
  if (fit$rank < ncol(design)) {
    j <- fit$pivot[1:fit$rank]
    fitted.values <- fit$coefficients[, j, drop = FALSE] %*%
      t(fit$design[, j, drop = FALSE])
  }
  else {
    fitted.values <- fit$coefficients %*% t(fit$design)
  }
  fitted.cpm <- 2^fitted.values
  fitted.count <- 1e-06 * t(t(fitted.cpm) * (lib.size + 1))
  fitted.logcount <- log2(fitted.count)
  w <- 1/f(fitted.logcount)^4
  dim(w) <- dim(fitted.logcount)

  ## If array weights are provided, the values for `w` would then be used for a
  ## second round of estimation, and scaled by these weights
  aw <- c()
  if (!is.null(w0)){
    ## Use all defaults from this function, providing the precision weights
    aw <- arrayWeights(
      object = cpm, design = design, weights = w, var.design = NULL,
      var.group = NULL, prior.n = 10, method = "auto", maxiter = 50,
      tol = 1e-5, trace = FALSE
    )
    w <- t(aw * t(w))
  }

  ## Initialise output
  nm <- colnames(cpm)
  if (is.null(nm)) nm <- seq_len(i)
  out <- list()
  out$E <- cpm
  out$weights <- w
  out$design <- design
  out$targets <- data.frame(
    sample = nm,
    lib.size = lib.size
  )
  if (!is.null(aw)) out$targets$sample.weights <- aw
  out$voom.xy <- list(
    x = sx, y = sy, Amean = fit$Amean,
    xlab = "log2( count size + 0.5 )", ylab = "Sqrt ( standard deviation )"
  )
  out$voom.line <- l

  ## Return an EList
  new("EList", out)

}
