#' @title Plot PCA For Any Assay
#'
#' @description Plot PCA for any assay within a SummarizedExperiment
#'
#' @details
#' Uses ggplot2 to create a PCA plot for the selected assay. Any numerical
#' transformation prior to performing the PCA can be specified using the
#' `trans` argument
#'
#' @return
#' A ggplot2 object
#'
#' @param x An object containing an assay slot
#' @param assay The assay to perform PCA on
#' @param colour The column name to be used for colours
#' @param shape The column name to be used for determining the shape of points
#' @param label The column name to be used for labels
#' @param show_points logical(1). Display the points. If `TRUE` any labels will
#' repel. If `FALSE`, labels will appear at the exact points
#' @param pc_x numeric(1) The PC to plot on the x-axis
#' @param pc_y numeric(1) The PC to plot on the y-axis
#' @param trans character(1). Any transformative function to be applied to the
#' data before performing the PCA, e.g. `trans = "log2"`
#' @param ... Passed to \link[stats]{prcomp}
#'
#' @examples
#' nrows <- 200; ncols <- 4
#' counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
#' df <- DataFrame(treat = c("A", "A", "B", "B"))
#' se <- SummarizedExperiment(
#'   assays = SimpleList(counts = counts),
#'   colData = df
#' )
#' plotAssayPCA(se, "counts", colour = "treat", label = "sample")
#'
#' @importFrom SummarizedExperiment colData assay
#' @importFrom tibble rownames_to_column
#' @importFrom broom tidy
#' @importFrom dplyr left_join filter rename
#' @importFrom tidyr pivot_wider
#' @importFrom scales percent
#' @importFrom stats prcomp
#' @importFrom ggplot2 ggplot aes_string geom_point geom_text labs
#' @importFrom ggrepel geom_text_repel
#'
#' @name plotAssayPCA
#' @rdname plotAssayPCA-methods
#' @export
#'
setGeneric(
  "plotAssayPCA",
  function(x, ...){standardGeneric("plotAssayPCA")}
)
#' @importFrom SummarizedExperiment assay
#'
#' @rdname plotAssayPCA-methods
#' @export
setMethod(
  "plotAssayPCA",
  signature = signature(x = "SummarizedExperiment"),
  function(
    x, assay = "counts",
    colour = c(), shape = c(), label = c(), show_points = TRUE,
    pc_x = 1, pc_y = 2, trans = c(), ...) {

    PC <- value <- c() # avoiding R CMD check errors

    if (is.null(colnames(x))) colnames(x) <- as.character(seq_len(ncol(x)))
    df <- as.data.frame(colData(x))
    df <- rownames_to_column(df, "sample")
    if (!is.null(colour)) stopifnot(colour %in% colnames(df))
    if (!is.null(shape)) stopifnot(shape %in% colnames(df))
    if (!is.null(label)) stopifnot(label %in% colnames(df))
    stopifnot(is.logical(show_points))
    stopifnot(length(pc_x) == 1 & length(pc_y) == 1)

    mat <- assay(x, assay)
    if (!is.null(trans)) mat <- match.fun(trans)(mat)

    pca <- prcomp(t(mat), ...)
    pca_df <- broom::tidy(pca)
    pca_df <- dplyr::rename(pca_df, sample = row)
    pca_df <- left_join(pca_df, df)
    pca_df <- dplyr::filter(pca_df, PC %in% c(pc_x, pc_y))
    pca_df <- pivot_wider(
      data = pca_df, names_from = PC, values_from = value, names_prefix = "PC"
    )
    prop_var <- pca$sdev^2/sum(pca$sdev^2)
    labs <- lapply(
      list(x = pc_x, y = pc_y),
      function(x) {
        paste0("PC ", x, "(", percent(prop_var[x], accuracy = 0.1), ")")
      }
    )
    col_x <- paste0("PC", pc_x)
    col_y <- paste0("PC", pc_y)

    p <- ggplot(
      pca_df, aes_string(col_x, col_y, colour = colour, shape = shape)
    )
    if (show_points) p <- p + geom_point()
    if (!is.null(label) & show_points)
      p <- p + geom_text_repel(aes_string(label = label), show.legend = FALSE)
    if (!is.null(label) & !show_points)
      p <- p + geom_text(aes_string(label = label), show.legend = FALSE)
    p + labs(
      x = labs$x, y = labs$y, shpae = shape, colour = colour
    )
  }
)
