#' @title Draw a heatmap from a single SummarizedExperiment assay
#'
#' @description Use ggplot2 to create a heatmap from a SummarizedExperiment object
#'
#' @details
#' Draw a heatmap containing selected values from an assay within a
#' SummarizedExperiment object. Columns within the colData element of the object
#' can be used to facet along the x-axis (e.g. treatment groups).
#' The maximum number of points is set to be 100, although this can be changed
#' easily should the plot require more ranges to be drawn.
#'
#' The averages across any grouping of samples can be drawn as a line plot on
#' the side of the y-axis by setting `ysideline = TRUE`, with groups as
#' specified in `yside_col`. This feature is added for the specific context of
#' neighbouring or overlapping ranges, and as such may be less informative in
#' any other scenario
#'
#' The returned object is a ggplot2 object so scales can easily be added after
#' heatmap creation using scale_fill_\* for the main heatmap, and
#' scale_colour_\* for any groupings along the y-axis
#'
#' @param x a SummarizedExperiment object
#' @param assay the assay to take values from
#' @param by_x the parameter to use for the x-axis. Will default to column names
#' but should be one value per sample, such as an additional column containing
#' shortened sample labels.
#' @param facet_x column from colData(x) which will be used to group samples
#' along the x-axis
#' @param ysideline logical(1) Draw a line across the side of the y-axis
#' summarising values for each range
#' @param yside_col column from colData(x) to group and colour the lines drawn
#' on the side of the y-axis. If grouping by treatment or replicate, the mean
#' values will be shown
#' @param trans character(1). Any transformative function to be applied to the
#' data before calculating the density, e.g. `trans = "log2"`
#' @param n_max Maximum number of ranges to draw
#' @param ... Not used
#'
#' @return
#' A `ggplot2` object. Scales and labels can be added using conventional
#' `ggplot2` syntax.
#'
#' @examples
#' nrows <- 10; ncols <- 4
#' counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
#' colnames(counts) <- paste0("Sample_", seq_len(ncols))
#' df <- DataFrame(treat = c("A", "A", "B", "B"))
#' se <- SummarizedExperiment(
#'   assays = SimpleList(counts = counts),
#'   colData = df
#' )
#' rowRanges(se) <- GRanges(paste0("chr1:", seq_len(nrows)))
#' plotAssayHeatmap(se, facet_x = "treat")
#'
#'
#' @name plotAssayHeatmap
#' @rdname plotAssayHeatmap-methods
#' @export
#'
setGeneric(
    "plotAssayHeatmap", function(x, ...) standardGeneric("plotAssayHeatmap")
)
#'
#' @import SummarizedExperiment
#' @importFrom tidyr unnest pivot_longer
#' @importFrom tidyselect everything all_of
#' @importFrom dplyr group_by summarise
#' @importFrom methods as
#' @importFrom rlang sym syms .data !!!
#' @importFrom forcats fct_inorder
#' @import ggside
#' @import ggplot2
#'
#' @rdname plotAssayHeatmap-methods
#' @export
setMethod(
    "plotAssayHeatmap",
    signature = signature(x = "SummarizedExperiment"),
    function(
        x, assay = "counts", by_x = "colnames", facet_x = NULL,
        ysideline = FALSE, yside_col = NULL, trans = NULL, n_max = 100, ...
    ) {

        ## Check column names & set plot aesthetics
        if (is.null(colnames(x))) colnames(x) <- as.character(seq_len(ncol(x)))
        col_data <- as_tibble(colData(x), rownames = "colnames")
        args <- colnames(col_data)
        msg <- "Any columns named 'colnames' or 'value' will be overwritten"
        if (any(colnames(colData(x)) %in% c("colnames", "value"))) message(msg)
        by_x <- match.arg(by_x, c("colnames", args))
        if (!is.null(facet_x)) facet_x <- match.arg(facet_x, args)
        if (ysideline & !is.null(yside_col))
            yside_col <- sym(match.arg(yside_col, c("colnames", args)))

        msg <- paste(
            "Only", n_max, "ranges can be drawn.",
            "Please change the n_max parameter if you wish to draw more."
        )
        if (nrow(x) > n_max) stop(msg)

        gr <- rowRanges(x)
        if (!is.null(gr)) {
            gr <- granges(gr)
            mcols(gr) <- assay(x, assay)
            df <- as_tibble(gr)
            by_y <- "range"
        } else {
            mat <- assay(x, assay)
            if (is.null(rownames(mat))) rownames(mat) <- seq_len(nrow(mat))
            df <- as_tibble(mat, rownames = "index")
            by_y <- "index"
        }
        df <- pivot_longer(
            df, cols = all_of(colnames(x)), names_to = "colnames"
        )
        df <- left_join(df, col_data, by = "colnames")
        ## Sort out any transformation
        if (!is.null(trans)) {
            fn <- match.fun(trans)
            test <- fn(1)
            trans_ok <- length(test) == 1
            if (!trans_ok) stop("This transformation is not applicable")
            df[["value"]] <- fn(df[["value"]])
        }
        fill_lab <- ifelse(is.null(trans), assay, paste(trans, assay))

        df[[by_x]] <- fct_inorder(df[[by_x]])
        df[[by_y]] <- fct_inorder(df[[by_y]])
        p <- ggplot(
            df, aes(x = !!sym(by_x), y = !!sym(by_y), fill = .data[["value"]])
        ) +
            geom_raster() +
            labs(fill = fill_lab) +
            scale_x_discrete(expand = expansion(0)) +
            scale_y_discrete(expand = expansion(0))
        if (ysideline) {
            syms <- c(as.character(by_y), as.character(yside_col))
            merged_df <- group_by(df, !!!syms(syms))
            merged_df <-
                summarise(merged_df, value = mean(!!sym("value")), .groups = "drop")
            p <- p +
                geom_ysideline(
                    aes(
                        x = !!sym("value"), y = as.integer(!!sym(by_y)),
                        colour = {{ yside_col}}, group = {{ yside_col }}
                    ), data = merged_df, orientation = "y"
                ) +
                ggside(collapse = "y")
        }
        if (!is.null(facet_x)) {
            fm <- as.formula(paste(".~", facet_x))
            p <- p + facet_grid(fm, scales = "free_x")
        }
        p

    }
)
