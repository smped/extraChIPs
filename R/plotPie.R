#' @title Draw Pie Graphs based on one or more columns
#'
#' @description Draw Pie Graphs based one or more data.frame columns
#'
#' @details
#' Using a `data.frame` as input, this function will draw pie graphs based
#' on one ore more columns, by simply counting the values in combination
#' across these columns.
#' One column must be selected for the fill as a bare minimum, with up to three
#' being possible.
#' Additional columns can be set for the x-axis to draw a series of pie-graphs
#' in a row, with a further column able to added to layout a series of pie
#' graphs in a grid
#'
#' If only one column/category is chosen, category labels will be added around
#' the edge of the plot
#'
#' If `show_total = TRUE` the overall counts for each pie graph will be added
#' in the centre using \link[ggplot2]{geom_label}.
#' Parameters for these labels are customisable
#'
#' @return
#' A ggplot2 object able to be customised with colour scales and themes.
#'
#' Also note that the $data element of the returned object will contain the
#' data.frame used for plotting. The additional column `label_radians`
#' represents the mid-point of each pie slice and can be used for manually
#' adding labels to each pie.
#' Only applies when plotting across the `x` or `y` axes
#'
#' @param object An object (`data.frame`)
#' @param fill The category/column used to fill the slices of the pie charts
#' @param x The second (optional) category/column to place along the x-axis
#' @param y The final (optional) category/column to plce along the y-axis
#' @param scale_by Scale the counts by this column. In this case of a GRanges
#' object this defaults to the count (scale_by = "n") but can also be specified
#' as being width of each range (scale_by = "width"). If choosing width, width
#' will be displayed in Kb
#' @param scale_factor When scaling by another column, such as width, totals
#' will be divided by this value, with 1000 being the default to provide output
#' in kb.
#' @param total_geom The geom_* to use for the totals at the centre of each pie.
#' Setting this to 'none' will disable totals
#' @param total_glue \link[glue]{glue} syntax to use for the totals in the
#' centre of each pie. The column 'N' will produce the totals and any other
#' values or formatting may be added here.
#' @param total_colour,total_fill,total_alpha,total_size Colour, fill, alpha and
#' size for the main totals in the centre of each pie chart
#' @param width Scale the width of all pies
#' @param min_p The minimum proportion of the total required for adding labels.
#' Effectively removes labels from pie charts with few members. Alternatively
#' when only one column is specified, categories below this will not be shown
#' around the edge of the plot
#' @param max_p only display labels for segments representing less than this
#' proportion of the total.
#' @param cat_geom The geom_* to use for category labels corresponding to each
#' slice of the pie. Setting this to 'none' will disable category labels
#' @param cat_glue \link[glue]{glue} syntax to use for the category labels
#' corresponding to each slice of the pie charts. The columns 'n' and 'p' can
#' be used to print totals and proportions for each slice.
#' @param cat_colour,cat_fill,cat_size,cat_alpha Colour, fill, size and alpha
#' for category labels
#' @param cat_adj Adjust category labels
#' @param hole_width Add a hole in the middle to turn the plot into a donut.
#' Values between zero and 1 work best. Only implemented for pie charts using
#' one value (i.e. fill)
#' @param ... Not used
#'
#' @examples
#' set.seed(200)
#' df <- data.frame(
#'   feature = sample(
#'     c("Promoter", "Enhancer", "Intergenic"), 200, replace = TRUE
#'   ),
#'   TF1 = sample(c("Up", "Down", "Unchanged"), 200, replace = TRUE),
#'   TF2 = sample(c("Up", "Down", "Unchanged"), 200, replace = TRUE),
#'   w = rexp(200)
#' )
#' plotPie(df, fill = "feature", total_glue = "N = {comma(N)}")
#' plotPie(
#'   df, fill = "feature", scale_by = "w", total_geom = "none",
#'   cat_glue = "{percent(p)}", cat_size = 5
#' )
#' plotPie(df, fill = "feature", x = "TF1")
#' plotPie(
#'   df, fill = "feature", x = "TF1", y = "TF2", min_p = 0.02,
#'   total_geom = "none", cat_glue = "{n} / {N}"
#'  ) +
#'  scale_fill_viridis_d() +
#'  theme_bw()
#'
#'
#' ## And using a GRanges object
#' data("ex_prom")
#' gr <- ex_prom
#' mcols(gr) <- df[seq_along(gr),]
#' ## Show values by counts
#' plotPie(gr, fill = "feature", total_size = 5)
#' ## Show values scaled by width of each range as a donut plot
#' plotPie(
#'   gr, fill = "feature", scale_by = "width", total_glue = "{round(N, 1)}kb",
#'   cat_glue = "{percent(p, 0.1)}", cat_size = 4, total_size = 5, hole_width = 0.2
#' )
#'
#' @importClassesFrom GenomicRanges GRanges
#' @rdname plotPie-methods
#' @export
setMethod(
    "plotPie",
    signature = signature(object = "GRanges"),
    function(object, scale_by = c("n", "width"), ...){
        df <- as.data.frame(object)
        scale_by <- match.arg(scale_by)
        if (scale_by == "n") df[["n"]] <- 1
        ## Set width to be in Kb by default
        if (scale_by == "width") df[["width"]] <- width(object)
        plotPie(df, scale_by = scale_by, ...)
    }
)
#' @importClassesFrom S4Vectors DataFrame
#' @rdname plotPie-methods
#' @export
setMethod(
    "plotPie",
    signature = signature(object = "DataFrame"),
    function(object, ...){
        object <- as.data.frame(object)
        plotPie(object, ...)
    }
)
#'
#' @rdname plotPie-methods
#' @importFrom rlang ensym
#' @export
setMethod(
    "plotPie",
    signature = signature(object = "data.frame"),
    function(
        object, fill, x, y, scale_by, scale_factor = 1e3, width = 0.8,
        total_geom = c("label", "text", "none"), total_glue = "{comma(N)}",
        total_colour = "black", total_fill = "white", total_alpha = 1,
        total_size = 3, min_p = 0.01, max_p = 1,
        cat_geom = c("label", "text", "none"),
        cat_glue = "{.data[[fill]]}\n{comma(n, 1)}\n({percent(p, 0.1)})",
        cat_colour = "black", cat_fill = "white", cat_size = 3, cat_alpha = 1,
        cat_adj = 0, hole_width = 0, ...
    ) {

        if (missing(fill)) stop("The initial category must be defined as 'fill'\n")
        total_geom <- match.arg(total_geom)
        cat_geom <- match.arg(cat_geom)
        fill <- as.character(ensym(fill))
        if (!missing(x)) x <- as.character(ensym(x))
        if (!missing(y)) y <- as.character(ensym(y))

        if (missing(x) & missing(y)) {
            p <- .plotSinglePie(
                df = object, fill = fill, width = width,
                .total_geom = total_geom, .total_glue = total_glue,
                .total_size = total_size, .total_colour = total_colour,
                .total_fill = total_fill,.total_alpha = total_alpha,
                .min_p = min_p, .max_p = max_p,
                .cat_geom = cat_geom, .cat_glue = cat_glue,
                .cat_colour = cat_colour, .cat_fill = cat_fill,
                .cat_size = cat_size, .cat_alpha = cat_alpha,
                .cat_adj = cat_adj,
                .scale_by = scale_by, .scale_factor = scale_factor,
                .hole_width = hole_width
            )
        }

        if (missing(x) & !missing(y))
            stop("Charts can only be drawn in rows when using two columns.\n")

        if (!missing(x) & missing(y))
            p <- .plotDoublePie(
                df = object, x = x, fill = fill, width = width,
                .total_geom = total_geom, .total_glue = total_glue,
                .total_size = total_size, .total_colour = total_colour,
                .total_fill = total_fill, .total_alpha = total_alpha,
                .min_p = min_p, .max_p = max_p, .cat_geom = cat_geom,
                .cat_glue = cat_glue, .cat_colour = cat_colour,
                .cat_fill = cat_fill, .cat_size = cat_size,
                .cat_alpha = cat_alpha, .cat_adj = cat_adj,
                .scale_by = scale_by, .scale_factor = scale_factor
            )

        if (!missing(x) & !missing(y))
            p <- .plotTriplePie(
                df = object, x = x, y = y, fill = fill, width = width,
                .total_geom = total_geom, .total_glue = total_glue,
                .total_size = total_size, .total_colour = total_colour,
                .total_fill = total_fill, .total_alpha = total_alpha,
                .min_p = min_p, .max_p = max_p, .cat_geom = cat_geom,
                .cat_glue = cat_glue, .cat_colour = cat_colour,
                .cat_fill = cat_fill, .cat_size = cat_size,
                .cat_alpha = cat_alpha, .cat_adj = cat_adj,
                .scale_by = scale_by, .scale_factor = scale_factor
            )

        p

    }
)

#' @importFrom dplyr group_by summarise mutate arrange
#' @importFrom glue glue
#' @importFrom scales comma percent
#' @importFrom rlang '!!' sym
#' @importFrom tidyr complete
#' @import ggplot2
.plotSinglePie <- function(
        df, fill, width, .total_geom, .total_glue, .total_size, .total_colour,
        .total_fill, .total_alpha, .min_p, .max_p, .cat_geom, .cat_glue,
        .cat_colour, .cat_fill, .cat_size, .cat_alpha, .text_width, .cat_adj,
        .scale_by, .scale_factor, .hole_width
) {

    stopifnot(fill %in% colnames(df))
    df[[fill]] <- as.factor(df[[fill]])
    x <- y <- lab <- c() # R CMD check error avoidance

    grp_df <- group_by(df, !!sym(fill))
    if (missing(.scale_by)) {
        summ_df <- summarise(grp_df, n = dplyr::n(), .groups = "drop")
    } else {
        stopifnot(is.numeric(df[[.scale_by]]) & is.numeric(.scale_factor))
        if (.scale_by == "n") .scale_factor <- 1
        summ_df <- summarise(grp_df, n = sum(!!sym(.scale_by)) / .scale_factor)
    }
    summ_df <- complete(summ_df, !!sym(fill), fill = list(n = 0))
    summ_df <- mutate(summ_df, p = n / sum(n), lab = glue(.cat_glue))
    summ_df <- arrange(summ_df, desc(!!sym(fill)))
    summ_df$y <- cumsum(summ_df$p)
    N <- sum(summ_df$n)

    x_lim <- c(-(width / 2 + .hole_width), width / 2 + .hole_width + .cat_adj)
    p <- ggplot(data = summ_df, aes(0, p, fill = !!sym(fill))) +
        geom_col(width = width) +
        xlim(x_lim)

    if (.cat_geom != "none") {
        geom <- paste0("geom_", .cat_geom)
        args <- list(
            mapping = aes(y = y - p / 2, label = lab),
            data = dplyr::filter(summ_df, p >= .min_p, p <= .max_p),
            x = .cat_adj,
            colour = .cat_colour, alpha = .cat_alpha, size = .cat_size,
            show.legend = FALSE
        )
        if (.cat_geom == "label") args$fill <- .cat_fill
        p <- p + do.call(geom, args)
    }

    if (.total_geom != "none") {
        lab_fun <- match.fun(paste0("geom_", .total_geom))
        p <- p + lab_fun(
            aes(x, y, label = lab),
            data = tibble(
                x = 0 - (width / 2 + .hole_width),
                y = 0, lab = glue(.total_glue)
            ),
            fill = .total_fill, colour = .total_colour,
            size = .total_size, alpha = .total_alpha,
            inherit.aes = FALSE
        )
    }

    p + coord_polar("y", start = 0) + theme_void()
}

#' @importFrom dplyr group_by summarise mutate ungroup distinct
#' @importFrom glue glue
#' @importFrom tidyr pivot_wider complete
#' @importFrom tidyselect all_of
#' @importFrom scales comma
#' @importFrom rlang '!!' sym .data
#' @importFrom ggforce stat_pie
#' @import ggplot2
.plotDoublePie <- function(
        df, x, fill, width, .total_geom, .total_glue, .total_size,
        .total_colour, .total_fill, .total_alpha, .min_p, .max_p, .cat_geom,
        .cat_glue, .cat_colour, .cat_fill, .cat_size, .cat_alpha, .cat_adj,
        .scale_by, .scale_factor
) {

    stopifnot(all(c(x, fill) %in% colnames(df)))
    df[[x]] <- as.factor(df[[x]])
    df[[fill]] <- as.factor(df[[fill]])
    r <- N <- value <- p <- c() # R CMD check error avoidance

    grp_df <- group_by(df, !!sym(x), !!sym(fill))
    if (missing(.scale_by)) {
        ## The column 'value' is hard-wired into geom_scatterpie for long format
        summ_df <- summarise(grp_df, value = dplyr::n(), .groups = "drop_last")
    } else {
        stopifnot(is.numeric(df[[.scale_by]]) & is.numeric(.scale_factor))
        if (.scale_by == "n") .scale_factor <- 1
        summ_df <- summarise(
            grp_df, value = sum(!!sym(.scale_by)) / .scale_factor,
            .groups = "drop_last"
        )
    }
    summ_df <- mutate(
        summ_df,
        p = value / sum(value), label_radians = 2 * pi * (cumsum(p) - 0.5 * p)
    )
    summ_df <- ungroup(mutate(summ_df, N = sum(value)))
    summ_df <- complete(
        summ_df, !!sym(fill), !!sym(x),
        fill = list(value = 0, N = 0, p = 0, label_radians = 0)
    )
    summ_df <- mutate(
        summ_df,
        n = value, lab = glue(.cat_glue),
        r = N / sum(N), r = 0.5 * r / max(r),
        x = as.integer(!!sym(x)),
        lab_x = x + 0.5 * r * sin(!!sym("label_radians")) * (1 + .cat_adj),
        lab_y = 1 + 0.5 * r * cos(!!sym("label_radians")) * (1 + .cat_adj)
    )

    p <- ggplot(data = summ_df) +
        ggforce::stat_pie(
            aes(
                x0 = x, y0 = 1, r0 = 0, r = width * r, fill = !!sym(fill),
                amount = value
            )
        ) +
        coord_equal() +
        scale_x_continuous(
            breaks = seq_along(levels(df[[x]])), labels = levels(df[[x]])
        ) +
        theme(
            panel.grid = element_blank(),
            axis.text.y = element_blank(), axis.title.y = element_blank(),
            axis.ticks.y = element_blank()
        )

    if (.total_geom != "none") {
        lab_fun <- match.fun(paste0("geom_", .total_geom))
        lab_df <- dplyr::filter(summ_df, N > .min_p * sum(N))
        lab_df <- distinct(lab_df, x, N, .keep_all = TRUE)
        lab_df <- mutate(lab_df, y = 1, lab = glue(.total_glue))
        p <- p + lab_fun(
            aes(x = .data[["x"]], y = .data[["y"]], label = .data[["lab"]]),
            data = lab_df,
            fill = .total_fill, colour = .total_colour, size = .total_size,
            alpha = .total_alpha, inherit.aes = FALSE
        )
    }

    if (.cat_geom != "none") {
        geom <- paste0("geom_", .cat_geom)
        args <- list(
            mapping = aes(
                .data[["lab_x"]], .data[["lab_y"]], label = .data[["lab"]]
            ),
            data = dplyr::filter(
                summ_df, n / sum(n) >= .min_p, n / sum(n) <= .max_p
            ),
            colour = .cat_colour, alpha = .cat_alpha, size = .cat_size,
            show.legend = FALSE
        )
        if (.cat_geom == "label") args$fill <- .cat_fill
        p <- p + do.call(geom, args)
    }

    p + labs(x = x, fill = fill)

}

#' @importFrom dplyr group_by summarise mutate ungroup distinct
#' @importFrom tidyr pivot_wider complete
#' @importFrom tidyselect all_of
#' @importFrom scales comma
#' @importFrom rlang '!!' sym
#' @importFrom ggforce stat_pie
#' @importFrom glue glue
#' @import ggplot2
.plotTriplePie <- function(
        df, x, y, fill, width, .total_geom, .total_glue, .total_size,
        .total_colour, .total_fill, .total_alpha, .min_p, .max_p, .cat_geom,
        .cat_glue, .cat_colour, .cat_fill, .cat_size, .cat_alpha, .cat_adj,
        .scale_by, .scale_factor
) {

    stopifnot(all(c(x, y, fill) %in% colnames(df)))
    df[[x]] <- as.factor(df[[x]])
    df[[y]] <- as.factor(df[[y]])
    df[[fill]] <- as.factor(df[[fill]])
    r <- N <- p <- value <- c() # R CMD check error avoidance

    grp_df <- group_by(df, !!sym(x), !!sym(y), !!sym(fill))
    if (missing(.scale_by)) {
        summ_df <- summarise(grp_df, value = dplyr::n(), .groups = "drop_last")
    } else {
        stopifnot(is.numeric(df[[.scale_by]]) & is.numeric(.scale_factor))
        if (.scale_by == "n") .scale_factor <- 1
        summ_df <- summarise(
            grp_df, value = sum(!!sym(.scale_by)) / .scale_factor,
            .groups = "drop_last"
        )
    }
    summ_df <- mutate(
        summ_df,
        p = value / sum(value), label_radians = 2 * pi * (cumsum(p) - 0.5 * p)
    )
    summ_df <- ungroup(mutate(summ_df, N = sum(value)))
    summ_df <- complete(
        summ_df, !!sym(fill), !!sym(x), !!sym(y),
        fill = list(value = 0, N = 0, p = 0, label_radians = 0)
    )
    summ_df <- mutate(
        summ_df,
        n = value, lab = glue(.cat_glue),
        r = N / sum(N), r = 0.5 * r / max(r),
        x = as.integer(!!sym(x)), y = as.integer(!!sym(y)),
        lab_x = x + 0.5 * r * sin(!!sym("label_radians")) * (1 + .cat_adj),
        lab_y = y + 0.5 * r * cos(!!sym("label_radians")) * (1 + .cat_adj)
    )

    p <- ggplot(data = summ_df) +
        ggforce::stat_pie(
            aes(
                x0 = x, y0 = y, r0 = 0, r = width * r, fill = !!sym(fill),
                amount = value
            )
        ) +
        coord_equal() +
        scale_x_continuous(
            breaks = seq_along(levels(df[[x]])), labels = levels(df[[x]])
        ) +
        scale_y_continuous(
            breaks = seq_along(levels(df[[y]])), labels = levels(df[[y]])
        )

    if (.total_geom != "none") {
        lab_fun <- match.fun(paste0("geom_", .total_geom))
        lab_df <- distinct(summ_df, x, y, N)
        lab_df <- dplyr::filter(lab_df, N > .min_p * sum(N))
        lab_df <- mutate(lab_df, lab = glue(.total_glue))
        p <- p + lab_fun(
            aes(!!sym(x), !!sym(y), label = !!sym("lab")),
            data = lab_df,
            fill = .total_fill, colour = .total_colour,
            size = .total_size, alpha = .total_alpha,
            inherit.aes = FALSE
        )
    }

    if (.cat_geom != "none") {
        geom <- paste0("geom_", .cat_geom)
        args <- list(
            mapping = aes(!!sym("lab_x"), !!sym("lab_y"), label = !!sym("lab")),
            data = dplyr::filter(
                summ_df, n / sum(n) >= .min_p, n / sum(n) <= .max_p
            ),
            colour = .cat_colour, alpha = .cat_alpha, size = .cat_size,
            show.legend = FALSE
        )
        if (.cat_geom == "label") args$fill <- .cat_fill
        p <- p + do.call(geom, args)
    }

    p + labs(x = x, y = y, fill = fill)

}
