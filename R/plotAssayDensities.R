#' @title Plot Densities For Any Assay
#'
#' @description Plot Densities for any assay within a SummarizedExperiment
#'
#' @details
#' Uses ggplot2 to create a density plot for all samples within the selected
#' assay
#'
#' @return
#' A ggplot2 object
#'
#' @examples
#' nrows <- 200; ncols <- 4
#' counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
#' df <- DataFrame(treat = c("A", "A", "B", "B"))
#' se <- SummarizedExperiment(
#'   assays = SimpleList(counts = counts),
#'   colData = df
#' )
#' plotAssayDensities(se, colour_by = "treat")
#'
#' @param x A SummarizedExperiment object
#' @param assay An assay within x
#' @param colour_by The column in colData to colour lines by
#' @param group_by Draw lines based on this column in colData. Defaults to
#' each sample
#' @param ... Passed to \link[stats]{density}
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
  function(x, assay = "counts", colour_by = c(), group_by, ...) {

    ## Checks
    stopifnot(assay %in% assayNames(x))
    if (is.null(colnames(x))) colnames(x) <- as.character(seq_len(ncol(x)))

    col_data <- as.data.frame(colData(x))
    col_data <- rownames_to_column(col_data, "sample")
    if (missing(group_by)) group_by <- "sample"
    stopifnot(group_by %in% colnames(col_data))
    if (!is.null(colour_by))
      stopifnot(colour_by %in% colnames(col_data))

    mat <- assay(x, assay)
    dens <- apply(mat, MARGIN = 2, density, ...)
    dens <- lapply(dens, function(d){list(tibble(x = d$x, y = d$y))})
    df <- as_tibble(dens)
    df <- pivot_longer(
      df,
      cols = everything(), names_to = "sample", values_to = "dens"
    )
    df <- left_join(df, col_data, by = "sample")
    df <- unnest(df, dens)

    ggplot(df, aes_string("x", "y", group = group_by, colour = colour_by)) +
      geom_line() +
      labs(
        x = assay, y = "Density", colour = colour_by
      )
  }
)
