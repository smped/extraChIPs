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
#' @param trans character(1). NUuerical transformation to apply to the data
#' prior to RLE calculation
#' @param ... Not used
#'
#' @examples
#' nrows <- 200; ncols <- 4
#' counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
#' df <- DataFrame(treat = c("A", "A", "B", "B"))
#' se <- SummarizedExperiment(
#'   assays = SimpleList(counts = counts),
#'   colData = df
#' )
#' plotAssayRle(se, "counts", fill = "treat")
#' plotAssayRle(se, "counts", fill = "treat", by_x = "treat")
#'
#' @importFrom SummarizedExperiment assay colData
#' @importFrom tidyr unnest
#' @importFrom tidyselect all_of
#' @importFrom dplyr group_by mutate ungroup
#' @importFrom rlang '!!' sym .data
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
        x, assay = "counts", colour = NULL, fill = NULL, rle_group = NULL,
        by_x = NULL, n_max = Inf, trans = NULL, ...
    ) {

        if (is.null(colnames(x))) colnames(x) <- as.character(seq_len(ncol(x)))
        df <- as.data.frame(colData(x))
        df$colnames <- colnames(x)
        args <- colnames(df)
        if (!is.null(colour)) colour <- sym(match.arg(colour, args))
        if (!is.null(fill)) fill <- sym(match.arg(fill, args))
        if (!is.null(rle_group)) rle_group <- sym(match.arg(rle_group, args))
        if (!is.null(by_x)) {
            by_x <- match.arg(by_x, args)
        } else {
            by_x <- "colnames"
        }

        n_max <- min(nrow(x), n_max)
        ind <- seq_len(n_max)
        if (n_max < nrow(x)) ind <- sample.int(nrow(x), n_max, replace = FALSE)

        mat <- assay(x, assay)[ind,]
        if (!is.null(trans)) {
            mat <- match.fun(trans)(mat)
            trans_ok <- all(
                is.matrix(mat), nrow(mat) == length(ind),
                colnames(mat) == colnames(x)
            )
            if (!trans_ok) stop("This transformation is not applicable")
        }

        df$vals <- split(t(mat), seq_len(ncol(x)))
        df <- unnest(df, all_of("vals"))
        if (!is.null(rle_group)) df <- group_by(df, !!rle_group)
        df <- mutate(df, rle = !!sym("vals") - median(!!sym("vals")))
        df <- ungroup(df)

        xlab <- ifelse(by_x == "colnames", "Sample", by_x)
        ggplot(
            df,
            aes(.data[[by_x]], .data[["rle"]], fill = {{ fill }}, colour = {{ colour }})
        ) +
            geom_boxplot() +
            labs(x = xlab, y = "RLE")

    }
)
