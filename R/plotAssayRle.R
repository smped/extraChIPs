#' @title Plot RLE for a given assay within a SummarizedExperiment
#'
#' @description Plot RLE for a given assay within a SummarizedExperiment
#'
#' @details
#' Uses ggplot2 to create an RLE plot for the selected assay. Any numerical
#' transformation prior to performing the RLE can be specified using the
#' `trans` argument
#'
#' @return
#' A ggplot2 object
#'
#' @param x A SummarizedExperiment object
#' @param assay The assay to plot
#' @param colour Column from `colData(x)` to outline the boxplots
#' @param fill Column from `colData(x)` to fill the boxplots
#' @param rle_group Column from `colData(x)` to calculate RLE within groups
#' Commonly an alternative sample label.
#' @param by_x Boxplots will be drawn by this grouping variable from
#' `colData(x)`. If not specified, the default values will be `colnames(x)`
#' @param n_max Maximum number of points to plot
#' @param trans character(1). Numerical transformation to apply to the data
#' prior to RLE calculation
#' @param ... Passed to \link[ggplot2]{geom_boxplot}
#'
#' @examples
#' data("se")
#' se$treatment <- c("E2", "E2", "E2", "E2DHT", "E2DHT", "E2DHT")
#' se$sample <- colnames(se)
#' ## A conventional RLE Plot using all samples
#' plotAssayRle(se, trans = "log1p", fill = "treatment")
#' ## Calculate RLE within groups
#' plotAssayRle(se, trans = "log1p", fill = "treatment", rle_group = "treatment")
#' # Or show groups combined
#' plotAssayRle(se, trans = "log1p", fill = "treatment", by_x = "treatment")
#'
#' @import SummarizedExperiment
#' @importFrom tidyr unnest
#' @importFrom tidyselect all_of
#' @importFrom dplyr group_by mutate ungroup
#' @importFrom rlang !! sym .data ensym
#' @importFrom stats median
#' @import ggplot2
#'
#' @rdname plotAssayRle-methods
#' @aliases plotAssayRle
#' @export
setMethod(
    "plotAssayRle",
    signature = signature(x = "SummarizedExperiment"),
    function(
        x, assay = "counts", colour, fill, rle_group, by_x, n_max = Inf,
        trans = NULL, ...
    ) {

        if (is.null(colnames(x))) colnames(x) <- as.character(seq_len(ncol(x)))
        df <- as.data.frame(colData(x))
        df$colnames <- colnames(x)
        args <- colnames(df)
        if (missing(colour)) {
            colour <- NULL
        } else {
            colour <- as.character(ensym(colour))
            colour <- sym(match.arg(colour, args))
        }
        if (missing(fill)) {
            fill <- NULL
        } else {
            fill <- as.character(ensym(fill))
            fill <- sym(match.arg(fill, args))
        }
        if (missing(rle_group)) {
            rle_group <- NULL
        } else {
            rle_group <- as.character(ensym(rle_group))
            rle_group <- sym(match.arg(rle_group, args))
        }
        if (missing(by_x)) {
            by_x <- sym("colnames")
            x_lab <- "Sample"
        } else {
            by_x <- as.character(ensym(by_x))
            by_x <- sym(match.arg(by_x, args))
            x_lab <- as.character(by_x)
        }

        n_max <- min(nrow(x), n_max)
        ind <- seq_len(n_max)
        if (n_max < nrow(x)) ind <- sample.int(nrow(x), n_max, replace = FALSE)

        mat <- assay(x, assay)[ind,]
        y_lab <- paste0("RLE (", assay, ")")
        if (!is.null(trans)) {
            mat <- match.fun(trans)(mat)
            trans_ok <- all(
                is.matrix(mat), nrow(mat) == length(ind),
                colnames(mat) == colnames(x)
            )
            if (!trans_ok) stop("This transformation is not applicable")
            y_lab <- paste0("RLE (", trans, " ", assay, ")")
        }

        df$vals <- split(t(mat), seq_len(ncol(x)))
        df <- unnest(df, all_of("vals"))
        if (!is.null(rle_group)) df <- group_by(df, !!rle_group)
        df <- mutate(df, rle = !!sym("vals") - median(!!sym("vals")))
        df <- ungroup(df)

        ggplot(
            df,
            aes(!!by_x, .data[["rle"]], fill = !!fill, colour = !!colour)
        ) +
            geom_boxplot(...) +
            labs(x = x_lab, y = y_lab)

    }
)
