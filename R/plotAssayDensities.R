#' @title Plot Densities for any assay within a SummarizedExperiment
#'
#' @description Plot Densities for any assay within a SummarizedExperiment
#'
#' @details
#' Uses ggplot2 to create a density plot for all samples within the selected
#' assay
#'
#' @return
#' A `ggplot2` object. Scales and labels can be added using conventional
#' `ggplot2` syntax. (See example)
#'
#' @examples
#' nrows <- 200; ncols <- 4
#' counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
#' df <- DataFrame(treat = c("A", "A", "B", "B"))
#' se <- SummarizedExperiment(
#'   assays = SimpleList(counts = counts),
#'   colData = df
#' )
#' plotAssayDensities(se, colour = "treat") +
#'   labs(colour = "Treat")
#'
#' @param x A SummarizedExperiment object
#' @param assay An assay within x
#' @param group Draw lines based on this column in colData. Defaults to
#' each sample
#' @param colour The column in colData to colour lines by
#' @param linetype Any column in colData used to determine linetype
#' @param trans character(1). Any transformative function to be applied to the
#' data before calculating the density, e.g. `trans = "log2"`
#' @param n_max Maximum number of points to use when calculating densities
#' @param ... Not used
#'
#' @name plotAssayDensities
#' @rdname plotAssayDensities-methods
#' @export
#'
setGeneric(
  "plotAssayDensities",
  function(x, ...){standardGeneric("plotAssayDensities")}
)
#' @importFrom SummarizedExperiment assayNames colData
#' @importFrom stats density
#' @importFrom tibble as_tibble tibble rownames_to_column
#' @importFrom tidyr pivot_longer unnest
#' @importFrom tidyselect everything
#' @importFrom ggplot2 ggplot aes_string geom_line labs
#' @importFrom methods as
#' @importFrom dplyr left_join
#'
#' @rdname plotAssayDensities-methods
#' @export
setMethod(
  "plotAssayDensities",
  signature = signature(x = "SummarizedExperiment"),
  function(
    x, assay = "counts",
    group, colour = NULL, linetype = NULL, trans = NULL, n_max = Inf, ...
  ) {

    ## Checks
    stopifnot(assay %in% assayNames(x))
    if (is.null(colnames(x))) colnames(x) <- as.character(seq_len(ncol(x)))

    col_data <- as.data.frame(colData(x))
    col_data <- rownames_to_column(col_data, "sample")
    if (missing(group)) group <- "sample"
    group <- match.arg(group, colnames(col_data))
    if (!is.null(colour)) colour <- match.arg(colour, colnames(col_data))
    if (!is.null(linetype)) linetype <- match.arg(linetype, colnames(col_data))

    n_max <- min(nrow(x), n_max)
    ind <- seq_len(n_max)
    if (n_max < nrow(x)) ind <- sample.int(nrow(x), n_max, replace = FALSE)

    mat <- assay(x[ind,], assay)
    if (!is.null(trans)) {
      mat <- match.fun(trans)(mat)
      trans_ok <- all(
        is.matrix(mat), nrow(mat) == length(ind), colnames(mat) == colnames(x)
      )
      if (!trans_ok) stop("This transformation is not applicable")
    }

    dens <- apply(mat, MARGIN = 2, density)
    dens <- lapply(dens, function(d){list(tibble(x = d$x, y = d$y))})
    df <- as_tibble(dens)
    df <- pivot_longer(
      df, cols = everything(), names_to = "sample", values_to = "dens"
    )
    df <- left_join(df, col_data, by = "sample")
    df <- unnest(df, dens)

    ggplot(
      df,
      aes_string("x", "y", group = group, colour = colour, linetype = linetype)
    ) +
      geom_line() +
      labs(
        x = ifelse(is.null(trans), assay, paste(trans, assay)),
        y = "Density", colour = colour, linetype = linetype
      )
  }
)
