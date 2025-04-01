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
#' plotAssayDensities(se, colour = "treatment", group = "treatment")
#' ## Use a data transformation
#' plotAssayDensities(se, trans = "log1p", colour = "treat")
#'
#' @param x A SummarizedExperiment object
#' @param assay An assay within x
#' @param group Used by \link[ggplot2]{geom_line}. Defaults to the column names
#' but treatment groups can also be specified to summarise within groups
#' @param colour,fill,alpha Optional column in colData to set the respective
#' aesthetics. Can also be any valid colour specification as a fixed value or
#' a fixed alpha value
#' @param linetype,linewidth Any optional column in colData used to determine
#' linetype or linewidth. Can also be fixed values
#' @param trans character(1). Any transformative function to be applied to the
#' data before calculating the density, e.g. `trans = "log2"`
#' @param n_max Maximum number of points to use when calculating densities
#' @param ... Passed to \link[ggplot2]{geom_density}
#'
#' @name plotAssayDensities
#' @rdname plotAssayDensities-methods
#' @export
#'
setGeneric(
    "plotAssayDensities", function(x, ...) standardGeneric("plotAssayDensities")
)
#' @import SummarizedExperiment
#' @importFrom rlang sym .data ensym
#' @import ggplot2
#'
#' @rdname plotAssayDensities-methods
#' @export
setMethod(
    "plotAssayDensities",
    signature = signature(x = "SummarizedExperiment"),
    function(
        x, assay = "counts", group = NULL, colour = NULL, fill = NULL,
        linetype = NULL, linewidth = NULL, alpha = NULL, trans = NULL,
        n_max = Inf, ...
    ) {

        ## Check column names
        if (is.null(colnames(x))) colnames(x) <- as.character(seq_len(ncol(x)))
        args <- colnames(colData(x))
        if (any(args %in% c("colnames", "vals"))) {
            msg <- "Any columns named 'colnames' or 'vals' will be overwritten"
            warning(msg)
            colData(x)$colnames <- NULL
            colData(x)$vals <- NULL
        }

        ## Now handle all of the possible options for plotting params
        param_list <- list(...)
        if (!is.null(colour)) {
            colour <- colour[[1]]
            if (.validColour(colour)) {
                param_list$colour <- colour
            } else {
                colour <- as.character(ensym(colour))
                colour <- sym(match.arg(colour, args))
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
        if (!is.null(linetype)) {
            linetype <- linetype[[1]]
            valid_lty <- c(
                scales::linetype_pal()(13), "blank", "solid", "dashed",
                "dotted", "dotdash", "longdash", "twodash"
            )
            if (is.numeric(linetype) | linetype %in% valid_lty) {
                param_list$linetype <- linetype
            } else {
                linetype <- as.character(ensym(linetype))
                linetype <- sym(match.arg(linetype, args))
            }
        }
        if (!is.null(linewidth)) {
            linewidth <- linewidth[[1]]
            if (is.numeric(linewidth)) {
                param_list$linewidth <- linewidth
            } else {
                linewidth <- as.character(ensym(linewidth))
                linewidth <- sym(match.arg(linewidth, args))
            }
        }
        if (!is.null(alpha)) {
            alpha <- alpha[[1]]
            if (is.numeric(alpha)) {
                param_list$alpha <- min(1, abs(alpha))
            } else {
                alpha <- as.character(ensym(alpha))
                alpha <- sym(match.arg(alpha, args))
            }
        }
        if (is.null(group)) {
            ## Use colnames as the default if not specified
            group <- sym("colnames")
        } else{
            group <- as.character(ensym(group))
            group <- sym(match.arg(group, args))
        }

        df <- .assayToDf(x, assay, n_max, trans)
        xlab <- ifelse(is.null(trans), assay, paste(trans, assay))
        ggplot(
            df,
            aes(
                x = .data[["vals"]], group = {{ group }},
                colour = {{ colour }}, fill = {{ fill }}, alpha = {{ alpha }},
                linetype = {{ linetype }}, linewidth = {{ linewidth }}
            )
        ) + do.call("geom_density", param_list) + labs(x = xlab, y = "Density")
    }
)

#' @keywords internal
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect everything
#' @importFrom dplyr left_join
.assayToDf <- function(x, assay, n_max, trans) {

    col_data <- as.data.frame(colData(x))
    col_data$colnames <- colnames(x)

    ## Subsample if required
    n_max <- min(nrow(x), n_max)
    ind <- seq_len(n_max)
    if (n_max < nrow(x)) ind <- sample.int(nrow(x), n_max, replace = FALSE)
    mat <- assay(x[ind, ], assay)
    if (!is.null(trans)) {
        mat <- match.fun(trans)(mat)
        trans_ok <- all(
            is.matrix(mat), nrow(mat) == length(ind),
            colnames(mat) == colnames(x)
        )
        if (!trans_ok) stop("This transformation is not applicable")
    }
    df <- as.data.frame(mat)
    df <- pivot_longer(df, everything(), names_to = "colnames", values_to = "vals")
    left_join(df, col_data, by = "colnames")

}
#' @keywords internal
.validColour <- function(x) {
    is_rgb <- grepl("^#[0-9A-F]+$", x) & nchar(x) %in% c(4, 7, 9)
    is_named <- x %in% grDevices::colours()
    is_rgb | is_named | is.numeric(x)
}
