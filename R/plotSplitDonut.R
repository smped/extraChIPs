#' @title Draw Two-Level Donut Charts
#'
#' @description Create Donut charts based on one or two columns in a data frame
#'
#' @details
#' Using a data.frame or GRanges object, this function enables creation of a
#' Pie/Donut chart with an inner and outer ring. The function itself is
#' extremely flexible allowing for separate colour palettes in the inner and
#' outer rings, as well as highly customisable labels.
#'
#' Sections can be exploded using a value from the inner ring or outer ring
#' separately, or in combination by setting `explode_query = "AND"`.
#' Exploded sections can be shifted by expanding the radius (`explode_r`), or
#' along the x/y co-ordinates using `explode_x/y`, allowing for detailed
#' placement of sections.
#'
#' If only the inner palette is specified, segments in the outer ring will be
#' assigned the same colours as the inner segments, but with increased
#' transparency. Only a single legend will be drawn in this scenario. If an
#' outer palette is specified, both colour palettes are completely distinct
#' and two distinct legends will be drawn. The placement of these legends,
#' along with the larger donut plot, can be manually specified by providing a
#' layout as defined in \link[patchwork]{plot_layout}. Names are not required
#' on this layout, but may be beneficial for code reproducibility.
#'
#' The inner label denoting the total can also be heavily customised using the
#' \link[glue]{glue} syntax to present the calculated value `N` along with any
#' additional text, such as 'kb' if scaling GenomicRanges by width. The same
#' approach can be taken for the inner and outer labels, where totals are
#' held in the value `n`, proportions are held in the value `p` and the values
#' corresponding to each segment can be accessed using `.data[[inner]]` or
#' `.data[[outer]]`. Column titles can be added using `{inner}`/`{outer}`.
#' Values from the inner segments can be added to the outer
#' labels using this strategy enabling a wide variety of labelling approaches
#' to be utilised.
#'
#' @return
#' A patchwork object consisting of both ggplot2 objects and legend grobs
#'
#' @param object A `GRanges` or `data.frame`-like object
#' @param inner Column name to create the inner ring
#' @param outer Column name to create the outer ring, subset by the inner ring
#' @param scale_by Column to scale values by. If provided, values in this column
#' will be summed, instead of simply counting entries. Any label in the centre
#' of the plot will also reflect this difference
#' @param scale_factor When scaling by another column, such as width, totals
#' will be divided by this value, with 1000 being the default to provide output
#' in kb.
#' @param r_centre The radius of the hole in the centre. Setting to zero will
#' create a Pie chart
#' @param r_inner,r_outer The radii of the inner/outer rings
#' @param total_size Label size total number of entries in the centre of the
#' plot.
#' @param total_colour Label colour for the summary total in the centre
#' @param total_glue \link[glue]{glue}-syntax for formatting the total which
#' appears in the centre of the plot. Internally, the value `N` will be
#' calculated and as such, this value should appear within this argument.
#' @param inner_glue,outer_glue \link[glue]{glue}-syntax for formatting labels
#' which appear on each inner/outer segment Internally, the values `n` and `p`
#' will be calculated as totals and proportions of the total. As such, these
#' values can appear within this argument, as well as the fields described in
#' the details
#' @param total_label,inner_label,outer_label Can take values 'text', 'label'
#' or 'none'. If setting one the first two values, the labelling function
#' `geom_*` will be called, otherwise no label will be drawn
#' @param label_size,inner_label_size,outer_label_size Size of all text labels
#' @param label_alpha,inner_label_alpha,outer_label_alpha transparency for
#' labels
#' @param label_colour,inner_label_colour,outer_label_colour Takes any colour
#' specification, with the additional option of 'palette'. In this special case,
#' the same palette as is used for each segment will be applied.
#' @param min_p,inner_min_p,outer_min_p only display labels for segments
#' representing greater than this proportion of the total. If inner/outer values
#' are specified, the values in `min_p` will be ignored for that layer
#' @param max_p,inner_max_p,outer_max_p only display labels for segments
#' representing less than this proportion of the total. If inner/outer values
#' are specified, the values in `max_p` will be ignored for that layer
#' @param inner_pattern,outer_pattern Regular expressions which are combined
#' with max_p and min_p values for accurately choosing labels
#' @param inner_rotate,outer_rotate logical(1). Rotate labels for inner or outer
#' rings. This will be ignored by when setting the geom as "label".
#' See \link[ggplot2]{geom_text}
#' @param explode_inner,explode_outer Regular expressions from either the inner
#' or outer ring for which segments will be 'exploded'
#' @param explode_query Setting to AND and specifying values for both the inner
#' and outer ring will require matches in both categories
#' @param explode_x,explode_y Numeric values for shifting exploded values
#' @param explode_r Radius expansion for exploded values
#' @param nudge_r,inner_nudge_r,outer_nudge_r Radius expansion for labels
#' @param expand Passed to \link[ggplot2]{expansion} for both x and y axes.
#' Can be helpful if labels are clipped by plot limits
#' @param inner_palette Colour palette for the inner ring
#' @param outer_palette Optional colour palette for the outer ring
#' @param inner_legend,outer_legend logical(1). Show legends for either layer
#' @param outer_p_by Scale the proportions for outer segments by the complete
#' dataset, or within each inner segment
#' @param layout Passed to \link[patchwork]{plot_layout}
#' @param ... Not used
#'
#' @examples
#' set.seed(200)
#' df <- data.frame(
#'   feature = sample(
#'     c("Promoter", "Enhancer", "Intergenic"), 200, replace = TRUE
#'   ),
#'   TF1 = sample(c("Up", "Down", "Unchanged"), 200, replace = TRUE),
#'   TF2 = sample(c("Up", "Down", "Unchanged"), 200, replace = TRUE)
#' )
#' ## The standard plot
#' plotSplitDonut(df, inner = "TF1", outer = "TF2", inner_legend = FALSE)
#'
#' ## Adding an exploded section along with an outer palette & customisation
#' plotSplitDonut(
#'   df, inner = "TF1", outer = "feature", total_label = "none",
#'   inner_label_alpha = 0.5, r_centre = 0,
#'   outer_glue = "{.data[[outer]]}\n(n = {n})", outer_label = "text",
#'   explode_inner = "Up", explode_outer = "Prom|Enh",
#'   explode_query = "AND", explode_r = 0.4, outer_rotate = TRUE,
#'   inner_palette = hcl.colors(3, "Spectral", rev = TRUE),
#'   outer_palette = hcl.colors(3, "Cividis")
#' )
#'
#' @importClassesFrom GenomicRanges GRanges
#' @importFrom rlang !! sym
#' @rdname plotSplitDonut-methods
#' @export
setMethod(
    "plotSplitDonut",
    signature = signature(object = "GRanges"),
    function(object, scale_by = c("n", "width"), ...){
        df <- as_tibble(object)
        scale_by <- match.arg(scale_by)
        if (scale_by == "n") df[["n"]] <- 1
        ## Set width to be in Kb by default
        if (scale_by == "width") df[["width"]] <- width(object)
        plotSplitDonut(df, scale_by = !!sym(scale_by), ...)
    }
)
#' @importClassesFrom S4Vectors DataFrame
#' @rdname plotSplitDonut-methods
#' @export
setMethod(
    "plotSplitDonut",
    signature = signature(object = "DataFrame"),
    function(object, ...){
        object <- as_tibble(object)
        plotSplitDonut(object, ...)
    }
)
#'
#'
#'
#' @importFrom dplyr group_by summarise mutate ungroup bind_rows case_when
#' @importFrom patchwork plot_layout area plot_spacer
#' @importFrom rlang !! !!! sym syms .data ensym
#' @importFrom scales comma percent
#' @importFrom ggforce geom_arc_bar
#' @importFrom stringr str_replace_all
#' @importFrom forcats fct_relabel
#' @importFrom glue glue
#'
#' @rdname plotSplitDonut-methods
#' @export
setMethod(
    "plotSplitDonut",
    signature = signature(object = "data.frame"),
    function(
        object, inner, outer, scale_by, scale_factor = 1e3,
        r_centre = 0.5, r_inner = 1, r_outer = 1,
        total_glue = "{comma(N)}",
        total_size = 5, total_colour = "black",
        inner_glue = "{inner} {.data[[inner]]}\n{percent(p,0.1)}",
        outer_glue = "{outer} {.data[[outer]]}\n{percent(p,0.1)}",
        total_label = c("label", "text", "none"),
        inner_label = c("label", "text", "none"),
        outer_label = c("label", "text", "none"),
        label_alpha = 1, inner_label_alpha = NULL, outer_label_alpha = NULL,
        label_size = 3, inner_label_size = NULL, outer_label_size = NULL,
        label_colour = "black", inner_label_colour = NULL,
        outer_label_colour = NULL,
        min_p = 0.05, inner_min_p = NULL, outer_min_p = NULL,
        max_p = 1, inner_max_p = NULL, outer_max_p = NULL,
        inner_pattern = ".", outer_pattern = ".",
        inner_rotate = FALSE, outer_rotate = FALSE,
        explode_inner = NULL, explode_outer = NULL,
        explode_query = c("AND", "OR"), explode_x = 0, explode_y = 0,
        explode_r = 0,
        nudge_r = 0.5, inner_nudge_r = NULL, outer_nudge_r = NULL,
        expand = 0.1, inner_palette = NULL, outer_palette = NULL,
        inner_legend = TRUE, outer_legend = TRUE,
        outer_p_by = c("all", "inner"),
        layout = c(main = area(1, 1, 12, 12), lg1 = area(2, 12), lg2 = area(11, 12)),
        ...
    ) {

        ## R CMD check declarations
        x0 <- y0 <- x1 <- yend <- explode <- ring <- c()

        if (missing(inner)) stop('argument "inner" is missing')
        if (missing(outer)) stop('argument "outer" is missing')
        inner <- as.character(ensym(inner))
        outer <- as.character(ensym(outer))
        if (missing(scale_by)) {
            scale_by <- NULL
        } else {
            scale_by <- as.character(ensym(scale_by))
        }
        stopifnot(all(c(inner, outer, scale_by) %in% colnames(object)))
        stopifnot(inner != outer)
        object[[inner]] <- as.factor(object[[inner]])
        object[[outer]] <- as.factor(object[[outer]])
        explode_query <- match.arg(explode_query)
        total_label <- match.arg(total_label)
        inner_label <- match.arg(inner_label)
        outer_label <- match.arg(outer_label)
        outer_p_by <- match.arg(outer_p_by)
        if (is.null(explode_inner)) explode_inner <- "^$"
        if (is.null(explode_outer)) explode_outer <- "^$"
        if (any(c(explode_inner, explode_outer) == "^$")) explode_query <- "OR"
        if (is.null(inner_min_p)) inner_min_p <- min_p
        if (is.null(outer_min_p)) outer_min_p <- min_p
        if (is.null(inner_max_p)) inner_max_p <- max_p
        if (is.null(outer_max_p)) outer_max_p <- max_p
        if (is.null(inner_label_alpha)) inner_label_alpha <- label_alpha
        if (is.null(outer_label_alpha)) outer_label_alpha <- label_alpha
        if (is.null(inner_label_size)) inner_label_size <- label_size
        if (is.null(outer_label_size)) outer_label_size <- label_size
        if (is.null(inner_label_colour)) inner_label_colour <- label_colour
        if (is.null(outer_label_colour)) outer_label_colour <- label_colour
        if (is.null(inner_nudge_r)) inner_nudge_r <- nudge_r
        if (is.null(outer_nudge_r)) outer_nudge_r <- nudge_r
        stopifnot(is.character(inner_pattern) & is.character(outer_pattern))

        summ_df <- group_by(object, !!!syms(c(inner, outer)))
        if (is.null(scale_by)) {
            summ_df <- summarise(summ_df, n = dplyr::n(), .groups = "drop_last")
            N <- nrow(object)
        } else {
            stopifnot(is.numeric(object[[scale_by]]) & is.numeric(scale_factor))
            if (scale_by == "n") scale_factor <- 1 # Will happen with GRanges
            summ_df <- summarise(
                summ_df, n = sum(!!sym(scale_by)) / scale_factor,
                .groups = "drop_last"
            )
            N <- sum(object[[scale_by]])
        }
        inner_df <- summarise(summ_df, n = sum(n), .groups = "drop")
        inner_df <- mutate(
            inner_df,
            ring = "inner", p = n / sum(n),
            x = r_centre, x1 = x + r_inner,
            y = cumsum(c(0, p[seq_len(nrow(inner_df) - 1)])),
            yend = cumsum(p),
            angle =  90 + (0.5 * p - cumsum(p)) * 360,
            angle = ifelse(
                abs(!!sym("angle")) > 90, !!sym("angle") + 180, !!sym("angle")
            )
        )
        if (!inner_rotate) inner_df$angle <- 0
        lev_inner <- levels(inner_df[[inner]])

        outer_df <- ungroup(summ_df)
        outer_df <- mutate(
            outer_df,
            ring = "outer", p = n / sum(n),
            x = r_inner + r_centre, x1 = x + r_outer,
            y = cumsum(c(0, p[seq_len(nrow(outer_df) - 1)])),
            yend = cumsum(p),
            angle =  90 + (0.5 * p - cumsum(p)) * 360,
            angle = ifelse(
                abs(!!sym("angle")) > 90, !!sym("angle") + 180, !!sym("angle")
            )
        )
        if (!outer_rotate) outer_df$angle <- 0
        if (outer_p_by == "inner")
            outer_df <- mutate(outer_df, p = p / sum(p), .by = !!sym(inner))
        lev_outer <- levels(outer_df[[outer]])
        lev_all <- unique(c(lev_inner, lev_outer))

        plot_df <- bind_rows(inner_df, outer_df)
        plot_df <- mutate(
            plot_df,
            start = y * 2 * pi, end = yend * 2 * pi, mid = 0.5 * (start + end),
            explode = case_when(
                explode_query == "AND" & grepl(explode_inner, !!sym(inner)) &
                    grepl(explode_outer, !!sym(outer)) ~ TRUE,
                explode_query == "OR" & (
                    grepl(explode_inner, !!sym(inner)) |
                        grepl(explode_outer, !!sym(outer))
                ) ~ TRUE,
                TRUE ~ FALSE
            ),
            x0 = ifelse(explode, explode_x, 0),
            y0 = ifelse(explode, explode_y, 0),
            x = case_when(
                explode_query == "OR" & grepl(explode_outer, !!sym(outer)) &
                    grepl(explode_inner, !!sym(inner)) ~ x + 2 * explode_r,
                explode ~ x + explode_r,
                TRUE ~ x
            ),
            x1 = case_when(
                explode_query == "OR" & grepl(explode_outer, !!sym(outer)) &
                    grepl(explode_inner, !!sym(inner)) ~ x1 + 2 * explode_r,
                explode ~ x1 + explode_r,
                TRUE ~ x1
            ),
            "colour" = fct_relabel(
                !!sym(inner), str_replace_all, "(.+)", paste(inner, "\\1")
            ),
            alpha = case_when(
                ring == "inner" ~ 1,
                ring == "outer" ~
                    as.numeric(!!sym(outer)) / (
                        max(as.numeric(!!sym(outer)), na.rm = TRUE) + 1
                    )
            ),
            "lab" = ifelse(ring == "inner", glue(inner_glue), glue(outer_glue))
        )

        ## Setup the default palette for the inner ring
        if (is.null(inner_palette))
            inner_palette <- hcl.colors(length(lev_inner), "Viridis")
        if (is.null(names(inner_palette)))
            inner_palette <- setNames(
                inner_palette[seq_along(lev_inner)], lev_inner
            )
        ## Make up a legend for the inner ring
        lg1 <- plot_spacer()
        if (inner_legend) {
            inner_leg <- ggplot(
                tibble("{inner}" := factor(lev_inner, levels = lev_inner)),
                aes(x = 1, y = .data[[inner]], fill = .data[[inner]])
            ) +
                geom_raster() +
                scale_fill_manual(values = inner_palette)
            inner_gr <- ggplot_gtable(ggplot_build(inner_leg))
            ind <- vapply(
                inner_gr$grobs, function(x) x$name == "guide-box", logical(1)
            )
            lg1 <- inner_gr$grobs[[which(ind)[[1]]]]
        }
        ## Rename for the actual plotting palette
        names(inner_palette) <- paste(inner, names(inner_palette))
        full_palette <- inner_palette

        ## Repeat for the outer ring if required
        lg2 <- plot_spacer()
        if (!is.null(outer_palette)) {
            stopifnot(length(outer_palette) >= length(lev_outer))
            if (is.null(names(outer_palette)))
                outer_palette <- setNames(
                    outer_palette[seq_along(lev_outer)], lev_outer
                )
            if (outer_legend) {
                outer_leg <- ggplot(
                    tibble("{outer}" := factor(lev_outer, levels = lev_outer)),
                    aes(x = 1, y = .data[[outer]], fill = .data[[outer]])
                ) +
                    geom_raster() +
                    scale_fill_manual(values = outer_palette)
                outer_gr <- ggplot_gtable(ggplot_build(outer_leg))
                ind <- vapply(
                    outer_gr$grobs, function(x) x$name == "guide-box",
                    logical(1)
                )
                lg2 <- outer_gr$grobs[[which(ind)[[1]]]]
            }
            names(outer_palette) <- paste(outer, names(outer_palette))
            full_palette <- c(inner_palette, outer_palette)
            ## Tidy up the main data.frame for the two palette style
            plot_df[["alpha"]] <- 1
            ind <- !is.na(plot_df[[outer]])
            plot_df[["colour"]] <- as.character(plot_df[["colour"]])
            plot_df[["colour"]][ind] <- paste(outer, plot_df[[outer]][ind])
            plot_df[["colour"]] <- factor(
                plot_df[["colour"]], levels = names(full_palette)
            )
        }

        ## Create the basic plot
        plt <- ggplot(plot_df) + geom_arc_bar(
            aes(
                x0 = x0, y0 = y0, r0 = x, r = x1,
                start = .data[["start"]], end = .data[["end"]],
                fill = .data[["colour"]], alpha = .data[["alpha"]]
            )
        )

        ## Add the centre label if required
        if (total_label != "none") {
            lab_fun <- match.fun(paste0("geom_", total_label))
            plt <- plt + lab_fun(
                aes(x, y, label = !!sym("lab")),
                data = tibble(x = 0, y = 0, "lab" = glue(total_glue)),
                size = total_size, inherit.aes = FALSE, colour = total_colour
            )
        }

        ## Add segment labels
        if (inner_label != "none") {
            plt <- .addLabel(
                plt = plt, df = plot_df, label_type = inner_label,
                label_colour = inner_label_colour, label_size = inner_label_size,
                label_alpha = inner_label_alpha, min_p = inner_min_p,
                max_p = inner_max_p, pattern = inner_pattern,
                nudge_r = inner_nudge_r, r = r_inner, .x = "x", .ring = "inner"
            )
        }
        if (outer_label != "none") {
            plt <- .addLabel(
                plt = plt, df = plot_df, label_type = outer_label,
                label_colour = outer_label_colour, label_size = outer_label_size,
                label_alpha = outer_label_alpha, min_p = outer_min_p,
                max_p = outer_max_p, pattern = outer_pattern,
                nudge_r = outer_nudge_r, r = 1, .x = "x1", .ring = "outer"
            )
        }

        ## Add all remaining elements & formatting
        plt <- plt +
            coord_equal() +
            theme_void() +
            scale_x_continuous(expand = expansion(expand)) +
            scale_y_continuous(expand = expansion(expand)) +
            scale_fill_manual(values = full_palette) +
            scale_colour_manual(values = full_palette) +
            scale_alpha_continuous(limits = c(0, 1)) +
            labs(fill = NULL) +
            guides(alpha = "none", fill = "none")

        plt + lg1 + lg2 + plot_layout(design = layout)

    }

)
#' @keywords internal
#' @importFrom rlang sym !!
.addLabel <- function(
        plt, df, label_type, label_colour, label_size, label_alpha, min_p, max_p,
        pattern, nudge_r, r, .x, .ring
) {

    ## Filter df for key parameters
    df <- dplyr::filter(
        df,
        !!sym("ring") == .ring, grepl(pattern, !!sym("lab")),
        !!sym("p") >= min_p, !!sym("p") <= max_p
    )
    ## Add labels
    lab_fun <- match.fun(paste0("geom_", label_type))
    if (label_colour != "palette") {
        plt <- plt + lab_fun(
            aes(
                x = sin(!!sym("mid")) * (!!sym(.x) + nudge_r * r) + !!sym("x0") ,
                y = cos(!!sym("mid")) * (!!sym(.x) + nudge_r * r) + !!sym("y0"),
                angle = !!sym("angle"), label = !!sym("lab")
            ),
            data = df,
            size = label_size, alpha = label_alpha, colour = label_colour
        )
    } else {
        plt <- plt + lab_fun(
            aes(
                x = sin(!!sym("mid")) * (!!sym(.x) + nudge_r * r) + !!sym("x0"),
                y = cos(!!sym("mid")) * (!!sym(.x) + nudge_r * r) + !!sym("y0"),
                colour = !!sym("colour"), angle = !!sym("angle"),
                label = !!sym("lab")
            ),
            data = df,
            size = label_size, alpha = label_alpha, show.legend = FALSE
        )
    }
    plt

}
