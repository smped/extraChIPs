#' @title Plot PCA For any assay within a SummarizedExperiment
#'
#' @description Plot PCA for any assay within a SummarizedExperiment object
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
#' @param n_max Subsample the data to this many points before performing PCA
#' @param ... Passed to \link[ggplot2]{geom_text}
#'
#' @examples
#' nrows <- 200; ncols <- 4
#' counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
#' df <- DataFrame(treat = c("A", "A", "B", "B"), sample = seq_len(4))
#' se <- SummarizedExperiment(
#'   assays = SimpleList(counts = counts),
#'   colData = df
#' )
#' plotAssayPCA(se, "counts", colour = "treat", label = "sample")
#' plotAssayPCA(
#'   se, "counts", colour = "treat", label = "sample",
#'   inherit.aes = FALSE, size = 5
#' )
#'
#'
#' @name plotAssayPCA
#' @rdname plotAssayPCA-methods
#' @export
#'
setGeneric(
    "plotAssayPCA",
    function(x, ...){standardGeneric("plotAssayPCA")}
)
#' @importFrom SummarizedExperiment colData assay
#' @importFrom broom tidy
#' @importFrom dplyr left_join
#' @importFrom tidyr pivot_wider
#' @importFrom scales percent
#' @importFrom stats prcomp
#' @importFrom ggrepel geom_text_repel
#' @import ggplot2
#'
#' @rdname plotAssayPCA-methods
#' @export
setMethod(
    "plotAssayPCA",
    signature = signature(x = "SummarizedExperiment"),
    function(
        x, assay = "counts", colour = NULL, shape = NULL, label = NULL,
        show_points = TRUE, pc_x = 1, pc_y = 2, trans = NULL, n_max = Inf, ...
    ) {

        if (is.null(colnames(x))) colnames(x) <- as.character(seq_len(ncol(x)))
        df <- as.data.frame(colData(x))
        args <- colnames(df)
        df$row <- rownames(df) ## To match tidy(pca) later
        if (!is.null(colour)) colour <- sym(match.arg(colour, args))
        if (!is.null(shape)) shape <- sym(match.arg(shape, args))
        if (!is.null(label)) label <- sym(match.arg(label, args))
        stopifnot(is.logical(show_points))

        n_max <- min(nrow(x), n_max)
        ind <- seq_len(n_max)
        if (n_max < nrow(x)) ind <- sample.int(nrow(x), n_max, replace = FALSE)

        mat <- assay(x[ind,], assay)
        if (!is.null(trans)) {
            mat <- match.fun(trans)(mat)
            trans_ok <- all(
                is.matrix(mat), nrow(mat) == length(ind),
                colnames(mat) == colnames(x)
            )
            if (!trans_ok) stop("This transformation is not applicable")
        }

        PC <- c() # avoiding R CMD check errors
        pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)
        max_comp <- length(pca$sdev)
        if (max(c(pc_x, pc_y)) > max_comp)
            stop("The highest available PC is ", max_comp)
        pca_df <- tidy(pca)
        pca_df <- left_join(pca_df, df, by = "row")
        pca_df <- dplyr::filter(pca_df, PC %in% c(pc_x, pc_y))
        pca_df <- pivot_wider(
            data = pca_df, names_from = "PC", values_from = "value",
            names_prefix = "PC"
        )
        prop_var <- pca$sdev^2 / sum(pca$sdev^2)
        labs <- lapply(
            c(x = pc_x[[1]], y = pc_y[[1]]),
            function(x) {
                paste0("PC", x, " (", percent(prop_var[x], accuracy = 0.1), ")")
            }
        )
        x <- sym(paste0("PC", pc_x))
        y <- sym(paste0("PC", pc_y))

        plot_aes <- aes(
            x = {{ x }}, y = {{ y }}, colour = {{ colour }}, shape = {{ shape }}
        )
        p <- ggplot(pca_df, plot_aes)
        if (show_points) p <- p + geom_point()
        if (!is.null(label)) {
            lab_fun <- ifelse(show_points, geom_text_repel, geom_text)
            if (show_points) formals(lab_fun)$show.legend <- FALSE
            p <- p +
                lab_fun(
                    aes(x = {{ x }}, y = {{ y }}, label = {{ label }}), ...
                )
        }

        p + labs(x = labs$x, y = labs$y, shape = shape, colour = colour)
    }
)
