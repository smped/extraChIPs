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
#' @param shape,size The column name(s) to be used for determining the shape
#' or size of points
#' @param label The column name to be used for labels
#' @param show_points logical(1). Display the points. If `TRUE` any labels will
#' repel. If `FALSE`, labels will appear at the exact points
#' @param pc_x numeric(1) The PC to plot on the x-axis
#' @param pc_y numeric(1) The PC to plot on the y-axis
#' @param trans character(1). Any transformative function to be applied to the
#' data before performing the PCA, e.g. `trans = "log2"`
#' @param n_max Subsample the data to this many points before performing PCA
#' @param tol Any rows with variance below this value will be excluded prior to
#' passing to \link[stats]{prcomp}. All rows are scaled and centred by default
#' @param rank Passed to \link[stats]{prcomp}
#' @param ... Passed to \link[ggplot2]{geom_text}
#'
#' @examples
#' data("se")
#' se$treatment <- c("E2", "E2", "E2", "E2DHT", "E2DHT", "E2DHT")
#' se$sample <- colnames(se)
#' plotAssayPCA(se, trans = "log1p", colour = "treatment", label = "sample")
#' plotAssayPCA(
#'   se, trans = "log1p", colour = "treatment", label = "sample",
#'   size = totals / 1e3
#' )
#' plotAssayPCA(
#'   se, trans = "log1p", colour = "treatment", label = "sample",
#'   show_points = FALSE
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
#' @import SummarizedExperiment
#' @importFrom broom tidy
#' @importFrom dplyr left_join
#' @importFrom tidyr pivot_wider
#' @importFrom scales percent
#' @importFrom stats prcomp
#' @importFrom ggrepel geom_text_repel
#' @importFrom matrixStats rowSds
#' @importFrom rlang sym ensym enexpr !!
#' @import ggplot2
#'
#' @rdname plotAssayPCA-methods
#' @export
setMethod(
    "plotAssayPCA",
    signature = signature(x = "SummarizedExperiment"),
    function(
        x, assay = "counts", colour, shape, size, label, show_points = TRUE,
        pc_x = 1, pc_y = 2, trans = NULL, n_max = Inf,
        tol = sqrt(.Machine$double.eps), rank = NULL, ...
    ) {

        if (is.null(colnames(x))) colnames(x) <- as.character(seq_len(ncol(x)))
        df <- as.data.frame(colData(x))
        args <- colnames(df)
        df$row <- rownames(df) ## To match tidy(pca) later
        if (missing(colour)) {
            colour <- NULL
        } else {
            colour <- as.character(ensym(colour))
            colour <- sym(match.arg(colour, args))
        }
        if (missing(shape)) {
            shape <- NULL
        } else {
            shape <- as.character(ensym(shape))
            shape <- sym(match.arg(shape, args))
        }
        if (missing(size)) {
            size <- NULL
        } else {
            ## This may be passed as a manpulation of data
            size <- enexpr(size)
            if (is.character(size)) size <- ensym(size)
        }
        if (missing(label)) {
            label <- NULL
        } else {
            label <- as.character(ensym(label))
            label <- sym(match.arg(label, args))
        }
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
        if (!is.null(tol)) {
            keep_rows <- rowSds(mat) >= tol
            if (sum(keep_rows) == 0)
                stop("Values are constant across all ranges")
            mat <- mat[keep_rows,]
        }

        PC <- c() # avoiding R CMD check errors
        pca <- prcomp(
            x = t(mat), center = TRUE, scale. = TRUE, tol = tol, rank. = rank
        )
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
            x = !!x, y = !!y, colour = !!colour, shape = !!shape, size = !!size
        )
        p <- ggplot(pca_df, plot_aes) + xlab(labs$x) + ylab(labs$y)
        if (show_points) p <- p + geom_point()
        if (!is.null(label)) {
            lab_fun <- ifelse(show_points, geom_text_repel, geom_text)
            if (show_points) formals(lab_fun)$show.legend <- FALSE
            p <- p +
                lab_fun(
                    aes(x = {{ x }}, y = {{ y }}, label = {{ label }}), ...
                )
        }
        p
    }
)

