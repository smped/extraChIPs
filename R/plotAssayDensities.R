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
#' `ggplot2` syntax.
#'
#' @examples
#' data("se")
#' se$treatment <- c("E2", "E2", "E2", "E2DHT", "E2DHT", "E2DHT")
#' ## Plot individual samples
#' plotAssayDensities(se, colour = "treatment")
#' ## Plot combined within treatment groups
#' plotAssayDensities(se, colour = "treatment", group = NULL)
#' ## Use a data transformation
#' plotAssayDensities(se, trans = "log1p", colour = "treatment")
#'
#' @param x A SummarizedExperiment object
#' @param assay An assay within x
#' @param colour The column in colData to colour lines by. To remove any
#' colours, set this argument to `NULL`
#' @param linetype Any optional column in colData used to determine linetype
#' @param group Used by \link[ggplot2]{geom_line}. Defaults to the sample names
#' but setting to NULL will over-write this and only groups specified by colour
#' or linetype will be drawn
#' @param trans character(1). Any transformative function to be applied to the
#' data before calculating the density, e.g. `trans = "log2"`
#' @param n_max Maximum number of points to use when calculating densities
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
#' @import SummarizedExperiment
#' @importFrom stats density
#' @importFrom tidyr unnest
#' @importFrom tidyselect everything all_of
#' @importFrom dplyr group_by summarise
#' @importFrom methods as
#' @importFrom rlang sym syms .data '!!!'
#' @import ggplot2
#'
#' @rdname plotAssayDensities-methods
#' @export
setMethod(
    "plotAssayDensities",
    signature = signature(x = "SummarizedExperiment"),
    function(
        x, assay = "counts", colour = NULL, linetype = NULL, group,
        trans = NULL, n_max = Inf, ...
    ) {

        ## Check column names & set plot aesthetics
        if (is.null(colnames(x))) colnames(x) <- as.character(seq_len(ncol(x)))
        col_data <- as.data.frame(colData(x))
        args <- colnames(col_data)
        msg <-
            "Any columns named 'colnames', 'vals' or 'dens' will be overwritten"
        if (any(args %in% c("colnames", "vals", "dens"))) message(msg)
        if (!is.null(colour)) colour <- sym(match.arg(colour, args))
        if (!is.null(linetype)) linetype <- sym(match.arg(linetype, args))
        ## By default, the plot should draw densities by sample.
        ## However,this should be able to be turned off by setting group to NULL
        if (missing(group)) {
            ## Use colnames as the default if not specified
            col_data$colnames <- colnames(x)
            group <- sym("colnames")
        } else{
            if (!is.null(group)) group <- sym(match.arg(group, args))
        }

        ## Subsample if required
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
        col_data$vals <- split(t(mat), seq_len(ncol(x)))
        df <- unnest(col_data, all_of("vals"))
        vars <- vapply(c(group, colour, linetype), as.character, character(1))
        df <- group_by(df, !!!syms(unique(vars)))
        df <- summarise(
            df, dens = list(as.data.frame(density(!!sym("vals"))[c("x", "y")])),
            .groups = "drop"
        )
        df <- unnest(df, all_of("dens"))
        xlab <- ifelse(is.null(trans), assay, paste(trans, assay))
        ggplot(
            df,
            aes(
                x = .data[["x"]], y = .data[["y"]], group = {{ group }},
                colour = {{ colour }}, linetype = {{ linetype }}
            )
        ) +
            geom_line() +
            labs(
                x = xlab, y = "Density", colour = colour, linetype = linetype
            )

    }
)
