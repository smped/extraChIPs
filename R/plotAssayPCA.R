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
#' @param colour,size The column names to be used for colours and point/label
#' size. Can be fixed values (e.g. size = 3) and can also be a manipulation of
#' a column, e.g. colour = log10(totals)
#' @param shape,fill The column name(s) to be used for determining the shape
#' or fill colour of plotted points. Can also be a fixed value
#' @param label The column name to be used for labels. Will default to the
#' column names of the SummarizedExperiment
#' @param show_points logical(1). Display the points. If `TRUE` any labels will
#' repel. If `FALSE`, labels will appear at the exact points
#' @param pc_x numeric(1) The PC to plot on the x-axis
#' @param pc_y numeric(1) The PC to plot on the y-axis
#' @param trans character(1). Any transformative function to be applied to the
#' data before performing the PCA, e.g. `trans = "log2"`
#' @param n_max Subsample the data to this many points before performing PCA
#' @param tol,rank Passed to \link[stats]{prcomp}
#' @param ... Passed to \link[ggplot2]{geom_text} and \link[ggplot2]{geom_point}
#'
#' @examples
#' data("se")
#' se$treatment <- c("E2", "E2", "E2", "E2DHT", "E2DHT", "E2DHT")
#' se$sample <- colnames(se)
#' plotAssayPCA(se, trans = "log1p", colour = "treatment", label = "sample")
#' plotAssayPCA(
#'   se, trans = "log1p", colour = "treatment", label = "sample",
#'   size = log10(totals), shape = 17
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
setGeneric("plotAssayPCA", function(x, ...) standardGeneric("plotAssayPCA"))
#' @import SummarizedExperiment
#' @import ggplot2
#' @importFrom tidyr pivot_wider
#' @importFrom rlang sym ensym enexpr !! enquo
#' @importFrom ggrepel geom_text_repel
#'
#' @rdname plotAssayPCA-methods
#' @export
setMethod(
    "plotAssayPCA",
    signature = signature(x = "SummarizedExperiment"),
    function(
        x, assay = "counts", colour = NULL, shape = NULL, size = NULL, fill = NULL,
        label = "colnames", show_points = TRUE, pc_x = 1, pc_y = 2, trans = NULL,
        n_max = Inf, tol = sqrt(.Machine$double.eps), rank = NULL, ...
    ) {


        if (is.null(colnames(x))) colnames(x) <- as.character(seq_len(ncol(x)))
        df <- as.data.frame(colData(x))
        df$colnames <- colnames(x)
        args <- colnames(df)

        ## PCs must be realistic
        max_comp <- ncol(x)
        if (max(c(pc_x, pc_y)) > max_comp)
            stop("The highest available PC is ", max_comp)

        param_list <- list(...)
        if (!is.null(shape)) {
            shape <- shape[[1]]
            if (is.numeric(shape) & shape %in% seq(0, 127)) {
                param_list$shape <- shape
            } else {
                shape <- as.character(ensym(shape))
                shape <- sym(match.arg(shape, args))
            }
        }
        if (!is.null(fill)) {
            fill <- fill[[1]]
            if (.validColour(fill)) {
                param_list$fill <- fill
            } else {
                fill <- as.character(ensym(fill))
                fill <- sym(match.arg(fill, args))
            }
        }
        ## This may be passed as a manipulation of data
        colour_quo <- enquo(colour)
        if (length(ls(environment(colour_quo))) == 0) {
            ## This will be a non formula
            if (!is.null(colour)) {
                colour <- colour[[1]]
                if (.validColour(colour)) {
                    param_list$colour <- colour
                } else {
                    colour <- as.character(ensym(colour))
                    colour <- sym(match.arg(colour, args))
                }
            }
        } else {
            colour <- colour_quo
        }
        ## This may also be passed as a manipulation of data
        size_quo <- enquo(size)
        if (length(ls(environment(size_quo))) == 0) {
            ## This will be a non formula
            if (!is.null(size)) {
                if (is.numeric(size)) {
                    param_list$size <- size
                } else {
                    size <- as.character(ensym(size))
                    size <- sym(match.arg(size, args))
                }
            }
        } else {
            size <- size_quo
        }
        if (!is.null(label)) {
            label <- as.character(ensym(label))
            label <- sym(match.arg(label, args))
        }
        stopifnot(is.logical(show_points))

        pca <- .runAssayPCA(x, assay, trans, n_max, tol, rank)
        pc_x <- paste0("PC", pc_x)
        pc_y <- paste0("PC", pc_y)
        pca_df <- data.frame( ## Similar to broom:tidy
            row = rep(rownames(pca$x), each = ncol(pca$x)),
            PC = rep(colnames(pca$x), times = nrow(pca$x)),
            value = as.numeric(t(pca$x))
        )
        pca_df <- subset(pca_df, pca_df$PC %in% c(pc_x, pc_y))
        pca_df <- cbind(pca_df, df[pca_df$row,])
        pca_df <- pivot_wider(
            data = pca_df, names_from = "PC", values_from = "value"
        )
        perc_var <- paste0(round(100 * pca$sdev^2 / sum(pca$sdev^2), 1), "%")
        names(perc_var) <- paste0("PC", seq_along(perc_var))
        labs <- lapply(
            c(x = pc_x[[1]], y = pc_y[[1]]),
            \(x) paste0(x, " (", perc_var[[x]], ")")
        )
        plot_aes <- aes(
            x = !!sym(pc_x), y = !!sym(pc_y), colour = {{ colour }},
            shape = {{ shape }}, size = {{ size }}, label = {{ label }},
            fill = {{ fill }}
        )
        if (is.numeric(shape)) plot_aes$shape <- NULL
        if (is.numeric(size)) plot_aes$size <- NULL
        p <- ggplot(pca_df, plot_aes) + labs(x = labs$x, y = labs$y)
        if (show_points) p <- p + do.call("geom_point", param_list)
        if (!is.null(label)) {
            lab_fun <- ifelse(show_points, "geom_text_repel", "geom_text")
            if (show_points) param_list$show.legend <- FALSE
            param_list$shape <- NULL
            param_list$fill <- NULL
            p <- p + do.call(lab_fun, param_list)
        }
        p
    }
)

#' @keywords internal
#' @importFrom stats prcomp
#' @importFrom matrixStats rowSds
.runAssayPCA <- function(x, assay, trans, n_max, tol, rank) {
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
        if (sum(keep_rows) == 0) stop("Values are constant across all ranges")
        mat <- mat[keep_rows,]
    }

    prcomp(
        x = t(mat), center = TRUE, scale. = TRUE, tol = tol, rank. = rank
    )
}
