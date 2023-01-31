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
#' `.data[[outer]]`. Values from the inner segments can be added to the outer
#' labels using this strategy enabling a wide variety of labellinf approaches
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
#' @param r_centre The radius of the hole in the centre. Setting to zero will
#' create a Pie chart
#' @param r_inner,r_outer The radii of the inner/outer rings
#' @param total_size Label size total number of entries in the centre of the
#' plot. Set to `NA` to hide the label itself
#' @param total_glue \link[glue]{glue}-syntax for formatting the total which
#' appears in the centre of the plot. Internally, the value `N` will be
#' calculated and as such, this value should appear within this argument.
#' @param inner_glue,outer_glue \link[glue]{glue}-syntax for formatting labels
#' which appear on each inner/outer segment Internally, the values `n` and `p`
#' will be calculated as totals and proportions of the total. As such, these
#' values can appear within this argument.
#' @param inner_label,outer_label Can take values 'text', 'label' or 'none'.
#' If setting one the first two values, the labelling function `geom_*` will be
#' called, otherwise no label will be drawn
#' @param label_size Size of all text labels
#' @param label_alpha transparency for labels in the inner ring only
#' @param min_p only display labels for segments representing greater than this
#' proportion of the total
#' @param explode_inner,explode_outer Regular expressions from either the inner
#' or outer ring for which segments will be 'exploded'
#' @param explode_query Setting to AND and specifying values for both the inner
#' and outer ring will require matches in both categories
#' @param explode_x,explode_y Numeric values for shifting exploded values
#' @param explode_r Radius expansion for exploded values
#' @param nudge_r Radius expansion for labels in the outer ring
#' @param expand Passed to \link[ggplot2]{expansion} for both x and y axes
#' @param inner_palette Colour palette for the inner ring
#' @param outer_palette Optional colour palette for the outer ring
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
#' plotSplitDonut(df, inner = "TF1", outer = "TF2")
#'
#' ## Adding an exploded section along with an outer palette & customisation
#' plotSplitDonut(
#'   df, inner = "TF1", outer = "feature", total_size = NA,
#'   label_alpha = 0.5, r_centre = 0,
#'   outer_glue = "{.data[[outer]]}\n(n = {n})", outer_label = "text",
#'   explode_inner = "Up", explode_outer = "Prom|Enh",
#'   explode_query = "AND", explode_r = 0.4,
#'   inner_palette = hcl.colors(3, "Spectral", rev = TRUE),
#'   outer_palette = hcl.colors(3, "Cividis")
#' )
#'
#' @importClassesFrom GenomicRanges GRanges
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
        if (scale_by == "width") df[["width"]] <- width(object) / 1e3
        plotSplitDonut(df, scale_by = scale_by, ...)
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
#' @importFrom rlang '!!' '!!!' sym syms .data
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
        object, inner, outer, scale_by = NULL,
        r_centre = 0.5, r_inner = 1, r_outer = 1,
        total_size = 5, total_glue = "{comma(N)}",
        inner_glue = "{inner} {.data[[inner]]}\n{percent(p,0.1)}",
        outer_glue = "{outer} {.data[[outer]]}\n{percent(p,0.1)}",
        inner_label = c("label", "text", "none"),
        outer_label = c("label", "text", "none"),
        label_alpha = 1, label_size = 3, min_p = 0.05,
        explode_inner = NULL, explode_outer = NULL,
        explode_query = c("AND", "OR"), explode_x = 0, explode_y = 0,
        explode_r = 0, nudge_r = 0.5,
        expand = 0.1, inner_palette = NULL, outer_palette = NULL,
        layout = c(main = area(1, 1, 6, 6), lg1 = area(2, 7), lg2 = area(4, 7)),
        ...
    ) {

        inner <- inner[[1]]
        outer <- outer[[1]]
        scale_by <- scale_by[[1]]
        stopifnot(all(c(inner, outer, scale_by) %in% colnames(object)))
        stopifnot(inner != outer)
        if (!is.null(scale_by)) stopifnot(is.numeric(object[[scale_by]]))
        object[[inner]] <- as.factor(object[[inner]])
        object[[outer]] <- as.factor(object[[outer]])
        explode_query <- match.arg(explode_query)
        inner_label <- match.arg(inner_label)
        outer_label <- match.arg(outer_label)
        if (is.null(explode_inner)) explode_inner <- "^$"
        if (is.null(explode_outer)) explode_outer <- "^$"
        if (any(c(explode_inner, explode_outer) == "^$")) explode_query <- "OR"

        summ_df <- group_by(object, !!!syms(c(inner, outer)))
        if (is.null(scale_by)) {
            summ_df <- summarise(summ_df, n = dplyr::n(), .groups = "drop_last")
            N <- nrow(object)
        } else {
            summ_df <- summarise(
                summ_df, n = sum(!!sym(scale_by)), .groups = "drop_last"
            )
            N <- sum(object[[scale_by]])
        }
        inner_df <- summarise(summ_df, n = sum(n), .groups = "drop")
        inner_df <- mutate(
            inner_df,
            ring = "inner", p = n / sum(n),
            x = r_centre, x1 = x + r_inner,
            y = cumsum(c(0, p[seq_len(nrow(inner_df) - 1)])), yend = cumsum(p)
        )
        lev_inner <- levels(inner_df[[inner]])

        outer_df <- ungroup(summ_df)
        outer_df <- mutate(
            outer_df,
            ring = "outer", p = n / sum(n),
            x = r_inner + r_centre, x1 = x + r_outer,
            y = cumsum(c(0, p[seq_len(nrow(outer_df) - 1)])), yend = cumsum(p)
        )
        lev_outer <- levels(outer_df[[outer]])
        lev_all <- unique(c(lev_inner, lev_outer))

        explode <- c() # R CMD check
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
            colour = fct_relabel(
                !!sym(inner), str_replace_all, "(.+)", paste(inner, "\\1")
            ),
            alpha = case_when(
                ring == "inner" ~ 1,
                ring == "outer" ~
                    as.numeric(!!sym(outer)) / (
                        max(as.numeric(!!sym(outer)), na.rm = TRUE) + 1
                    )
            ),
            lab = ifelse(ring == "inner", glue(inner_glue), glue(outer_glue))
        )

        ## Setup the default palette for the inner ring
        if (is.null(inner_palette))
            inner_palette <- hcl.colors(length(lev_inner), "Viridis")
        if (is.null(names(inner_palette)))
            inner_palette <- setNames(
                inner_palette[seq_along(lev_inner)], lev_inner
            )
        ## Make up a legend for the inner ring
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
        lg2 <- plot_spacer()
        ## Rename for the actual plotting palette
        names(inner_palette) <- paste(inner, names(inner_palette))
        full_palette <- inner_palette

        ## Repeat for the outer ring if required
        if (!is.null(outer_palette)) {
            stopifnot(length(outer_palette) >= length(lev_outer))
            if (is.null(names(outer_palette)))
                outer_palette <- setNames(
                    outer_palette[seq_along(lev_outer)], lev_outer
                )
            outer_leg <- ggplot(
                tibble("{outer}" := factor(lev_outer, levels = lev_outer)),
                aes(x = 1, y = .data[[outer]], fill = .data[[outer]])
            ) +
                geom_raster() +
                scale_fill_manual(values = outer_palette)
            outer_gr <- ggplot_gtable(ggplot_build(outer_leg))
            ind <- vapply(
                outer_gr$grobs, function(x) x$name == "guide-box", logical(1)
            )
            lg2 <- outer_gr$grobs[[which(ind)[[1]]]]
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
        x <- y <- x0 <- y0 <- x1 <- lab <- ring <- mid <- yend <- colour <- c()
        p <- ggplot(plot_df) + geom_arc_bar(
            aes(
                x0 = x0, y0 = y0, r0 = x, r = x1, start = .data[["start"]],
                end = .data[["end"]], fill = colour, alpha = .data[["alpha"]]
            )
        )

        ## Add the centre total if required
        if (!is.na(total_size)) p <- p + geom_label(
            x = 0 , y = 0, size = total_size, label = glue(total_glue)
        )

        ## Add segment labels
        if (inner_label != "none") {
            lab_fun <- match.fun(paste0("geom_", inner_label))
            p <- p + lab_fun(
                aes(
                    x = sin(mid) * (x + 0.5 * r_inner) + x0,
                    y = cos(mid) * (x + 0.5 * r_inner) + y0, label = lab
                ),
                data = dplyr::filter(plot_df, ring == "inner", p > min_p),
                size = label_size, alpha = label_alpha
            )
        }
        if (outer_label != "none") {
            lab_fun <- match.fun(paste0("geom_", outer_label))
            p <- p + lab_fun(
                aes(
                    x = sin(mid) * (x1 + nudge_r) + x0 ,
                    y = cos(mid) * (x1 + nudge_r) + y0, label = lab
                ),
                data = dplyr::filter(plot_df, ring == "outer", p > min_p),
                size = label_size, alpha = label_alpha
            )
        }

        ## Add all remaining elements & formatting
        p <- p +
            coord_equal() +
            theme_void() +
            scale_x_continuous(expand = expansion(expand)) +
            scale_y_continuous(expand = expansion(expand)) +
            scale_fill_manual(values = full_palette) +
            scale_alpha_continuous(limits = c(0, 1)) +
            labs(fill = NULL) +
            guides(alpha = "none", fill = "none")

        p + lg1 + lg2 +  plot_layout(design = layout)

    }

)
