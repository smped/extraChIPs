#' @title Draw Two-Level Donut Charts
#'
#' @description Create Donut charts based on one or two columns in a data frame
#'
#' @details
#' Using a data.frame or GRanges object, this function enables creation of a
#' Pie/Donut chart with an inner and outer ring
#'
#' @return
#' A patchwork object consisting of both ggplot2 objects and legend grobs
#'
#' @param object A `GRanges` or `data.frame`-like object
#' @param inner Column name to create the inner ring
#' @param outer Column name to create the outer ring, subset by the inner ring
#' @param show_total Display the total number of entries in the centre of the
#' plot
#' @param r_centre The radius of the hole in the centre. Setting to zero will
#' create a Pie chart
#' @param r_inner,r_outer The radii of the inner/outer rings
#' @param show_cat logical(1) Show category names for each segment
#' @param show_n logical(1) Show category totals for each segment
#' @param show_percent logical(1) Show pecentages represented by each segment
#' @param label_width Width at which labels will try to wrap onto a new line
#' @param label_size Size of all text labels
#' @param sep character string to be used when concatenating categories, totals
#' and percentages
#' @param label_alpha transparency for labels in the inner ring only
#' @param min_p only display labels for segments representing greater than this
#' proportion of the total
#' @param explode_cat Single values from either the inner or outer ring for
#' which segments will be 'exploded'
#' @param explode_x,explode_y Numeric values for shifting exploded values
#' @param explode_r Radius expansion for exploded values
#' @param nudge_r Radius expansion for labels in the outer ring
#' @param expand Passed to \link[ggplot2]{expansion} for both x and y axes
#' @param inner_palette Colour palette for the inner ring
#' @param outer_palette Optional colour paletter for the outer ring
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
#'   df, inner = "TF1", outer = "feature",
#'   label_width = 1, show_n = FALSE, label_alpha = 0.5, r_centre = 0,
#'   explode_cat = "Up", explode_r = 0.2, explode_x = -0.2, explode_y = 0.2,
#'   inner_palette = hcl.colors(3, "Cividis"),
#'   outer_palette = hcl.colors(3, "Spectral")
#' )
#'
#' @importClassesFrom S4Vectors DataFrame
#' @rdname plotSplitDonut-methods
#' @export
setMethod(
    "plotSplitDonut",
    signature = signature(object = "DataFrame"),
    function(object, ...){
        object <- as.data.frame(object)
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
#' @importFrom stringr str_wrap
#' @importFrom ggforce geom_arc_bar
#'
#' @rdname plotSplitDonut-methods
#' @export
setMethod(
    "plotSplitDonut",
    signature = signature(object = "data.frame"),
    function(
        object, inner, outer, r_centre = 0.5, r_inner = 1, r_outer = 1,
        show_total = TRUE, show_cat = TRUE, show_n = TRUE, show_percent = TRUE,
        label_width = 10, sep = " ", label_alpha = 1, label_size = 4, min_p = 0.05,
        explode_cat = "", explode_x = 0, explode_y = 0, explode_r = 0,
        nudge_r = 0.5, expand = 0.1, inner_palette = NULL, outer_palette = NULL,
        layout = c(area(1, 1, 6, 6), area(2, 7), area(3, 7)),
        ...
    ) {

        stopifnot(all(c(inner, outer) %in% colnames(object)))
        object[[inner]] <- as.factor(object[[inner]])
        object[[outer]] <- as.factor(object[[outer]])

        summ_df <- group_by(object, !!!syms(c(inner, outer)))
        summ_df <- summarise(summ_df, n = dplyr::n(), .groups = "drop_last")
        inner_df <- summarise(summ_df, n = sum(n), .groups = "drop")
        inner_df <- mutate(
            inner_df,
            p = n / sum(n),
            x = r_centre, x1 = x + r_inner,
            y = cumsum(c(0, p[seq_len(nrow(inner_df) - 1)])), yend = cumsum(p),
            fill = !!sym(inner), ring = "inner"
        )
        lev_inner <- levels(inner_df[[inner]])

        outer_df <- ungroup(summ_df)
        outer_df <- mutate(
            outer_df,
            p = n / sum(n),
            x = r_inner + r_centre, x1 = x + r_outer,
            y = cumsum(c(0, p[seq_len(nrow(outer_df) - 1)])), yend = cumsum(p),
            fill = !!sym(outer), ring = "outer"
        )
        lev_outer <- levels(outer_df[[outer]])
        lev_all <- unique(c(lev_inner, lev_outer))

        fill <- c()
        plot_df <- bind_rows(inner_df, outer_df)
        plot_df <- mutate(
            plot_df,
            fill = factor(fill, levels = lev_all),
            start = y * 2 * pi, end = yend * 2 * pi, mid = 0.5 * (start + end),
            lab = case_when(
                show_cat & show_n & show_percent ~
                    paste(fill, comma(n, 1), percent(p), sep = sep),
                show_cat & show_n & !show_percent ~
                    paste(fill, comma(n, 1), sep = sep),
                show_cat & !show_n & show_percent ~
                    paste(fill, percent(p), sep = sep),
                !show_cat & show_n & show_percent ~
                    paste(comma(n, 1), percent(p), sep = sep),
                !show_cat & show_n & !show_percent ~ comma(n, 1),
                !show_cat & !show_n & show_percent ~ percent(p),
                show_cat & !show_n & !show_percent ~ as.character(fill),
                !show_cat & !show_n & !show_percent ~ ""
            ),
            lab = str_wrap(lab, width = label_width),
            x0 = case_when(
                !!sym(inner) == explode_cat | !!sym(outer) == explode_cat ~
                    explode_x,
                TRUE ~ 0
            ),
            y0 = case_when(
                !!sym(inner) == explode_cat | !!sym(outer) == explode_cat ~
                    explode_y,
                TRUE ~ 0
            ),
            x = case_when(
                (!!sym(inner) == explode_cat | !!sym(outer) == explode_cat) &
                    ring == "outer" ~ x + explode_r,
                TRUE ~ x
            ),
            x1 = case_when(
                !!sym(inner) == explode_cat | !!sym(outer) == explode_cat ~
                    x1 + explode_r,
                TRUE ~ x1
            ),
            colour = !!sym(inner),
            alpha = case_when(
                ring == "inner" ~ 1,
                ring == "outer" ~
                    as.numeric(!!sym(outer)) / (
                        max(as.numeric(!!sym(outer)), na.rm = TRUE) + 1
                    )
            )
        )

        ## Setup the default palette for the inner ring
        if (is.null(inner_palette))
            inner_palette <- hcl.colors(length(lev_inner), "Viridis")
        if (is.null(names(inner_palette)))
            inner_palette <- setNames(inner_palette[seq_along(lev_inner)], lev_inner)
        full_palette <- inner_palette
        inner_leg <- ggplot(
            tibble("{inner}" := factor(lev_inner, levels = lev_inner)),
            aes(x = 1, y = .data[[inner]], fill = .data[[inner]])
        ) +
            geom_raster() +
            scale_fill_manual(values = inner_palette)
        inner_gr <- ggplot_gtable(ggplot_build(inner_leg))
        ind <- vapply(inner_gr$grobs, function(x) x$name == "guide-box", logical(1))
        lg1 <- inner_gr$grobs[[which(ind)[[1]]]]
        lg2 <- plot_spacer()

        ## Add the palette for the outer ring if required
        if (!is.null(outer_palette)) {
            stopifnot(length(outer_palette) >= length(lev_outer))
            if (is.null(names(outer_palette)))
                outer_palette <- setNames(
                    outer_palette[seq_along(lev_outer)], lev_outer
                )
            full_palette <- c(inner_palette, outer_palette)[lev_all]
            plot_df[["alpha"]] <- 1
            plot_df[["colour"]] <- plot_df[["fill"]]
            outer_leg <- ggplot(
                tibble("{outer}" := factor(lev_outer, levels = lev_outer)),
                aes(x = 1, y = .data[[outer]], fill = .data[[outer]])
            ) +
                geom_raster() +
                scale_fill_manual(values = outer_palette)
            outer_gr <- ggplot_gtable(ggplot_build(outer_leg))
            ind <- vapply(outer_gr$grobs, function(x) x$name == "guide-box", logical(1))
            lg2 <- outer_gr$grobs[[which(ind)[[1]]]]
        }

        ## Create the basic plot
        x <- y <- x0 <- y0 <- x1 <- lab <- ring <- mid <- yend <- colour <- c()
        p <- ggplot(plot_df) +
            geom_arc_bar(
                aes(
                    x0 = x0, y0 = y0, r0 = x, r = x1,
                    start = .data[["start"]], end = .data[["end"]],
                    fill = colour, alpha = .data[["alpha"]]
                )
            )

        ## Add the centre total if required
        if (show_total)
            p <- p + geom_label(
                x = 0 , y = 0, size = label_size * 1.1,
                label = comma(nrow(object), 1)
            )

        ## Add segment labels
        if (any(c(show_cat, show_n, show_percent)))
            p <- p + geom_label(
                aes(
                    x = sin(mid) * (x + 0.5 * r_inner) + x0,
                    y = cos(mid) * (x + 0.5 * r_inner) + y0,
                    label = lab
                ),
                data = dplyr::filter(plot_df, ring == "inner", p > min_p),
                size = label_size, alpha = label_alpha
            ) +
            geom_text(
                aes(
                    x = sin(mid) * (x1 + nudge_r) + x0 ,
                    y = cos(mid) * (x1 + nudge_r) + y0, label = lab
                ),
                data = dplyr::filter(plot_df, ring == "outer", p > min_p),
                size = label_size
            )

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
