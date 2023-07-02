#' @title Plot Pairwise Values from a GRangeList
#'
#' @description
#' Plot Pairwise Values from a GRangeList by overlapping GRanges
#'
#' @details
#' This function enables pairwise plotting of two elements within a GRangesList.
#' All elements of the GRangesList will contain the same columns, so a set of
#' consensus ranges are first formed, before then taking all values from each
#' GRangesList element which overlap the range and producing a piarwise plot.
#'
#' Given that not all ranges will have values in both elements, side panels are
#' produced which can show the distribution of plotted values, along with those
#' which are only found in one of the foundational GRanges. These can take the
#' form of density, violin or boxplots.
#'
#' Addition columns, such as Differential Signal status can also be used to
#' form pairwise groups and colour the points.
#'
#' If a column in the GRangesList is suitable for labelling points, such as a
#' column with genes mapped to each range, this can be specified using the
#' argument `label = "col_to_label"`.
#' Only the furthest point from the origin will be labelled within each group
#' used to colour the points.
#' Labels will only be added if they lie beyond a circle of radius `min_d` from
#' the origin.
#' If multiple genes are mapped to the range, these will be separated by the
#' string provided in the `label_sep` argument.
#'
#' A regression line and correlation co-efficient are added to the plot by
#' default, but can be hidden easily if preferred
#'
#' @param x A GRangesList
#' @param var The colunm to compare between list elements
#' @param index Which list elements to compare
#' @param colour Optional column to use for combining across elements and
#' setting point colour
#' @param label Optional column to use for labelling ranges with the most
#' extreme changes
#' @param name_sep Character string to separate names of the GRangesList and
#' the selected column. Will appear as axis-labels
#' @param xside,yside Will call geom_(x/y)side* from the package ggside and
#' show additional panels on the top and right of the plot respectively
#' @param side_panel_width Set the relative widths of the side panels
#' @param side_alpha Set the transparency of any side_(x/y) panels
#' @param xside_axis_pos Position for axis_labels in the top panel when using
#' a discrete axis
#' @param yside_axis_label Wrapping for axis labels on the right-side panel when
#' using a discrete axis. Set to waiver() to turn off wrapping
#' @param line_col,line_type,line_width Parameters for adding a regression line
#' through the points. Set linecol to either `NULL` or `NA` to hide this line
#' @param rho_geom Used to add correlation coefficients for the two values
#' @param rho_pos Place the correlation coefficient within the plotting region
#' @param rho_col,rho_size,rho_alpha Parameters for displaying the correlation
#' @param label_geom Used to add labels from the column specified in label
#' @param label_width Label text will be truncated to this length
#' @param label_sep If multiple values (e.g. genes) are mapped to a range,
#' separate values using this string
#' @param label_size,label_alpha Passed to the geom used for adding labels
#' @param min_d Labels will only be added if the points lie circle beyond a
#' sircle of this radius
#' @param group_sep Text separator used to separate categories when specifying colour
#' @param simplify_equal logical(1) When combining columns from both elements
#' for the colour categories, should shared values be annotated as 'Both ...'
#' instead of having longer, more difficult to read annotations.
#' @param plot_theme Sets the initial theme by using the default theme for the
#' current R session via get_theme()
#' @param ... Passed to `geom_point()` for the main panel
#'
#' @return A `ggside` or `ggplot2` object
#'
#' @examples
#' theme_set(theme_bw())
#' set.seed(100)
#' gr1 <- GRanges(paste0("chr1:", seq(10, 150, by = 10)))
#' width(gr1) <- 5
#' gr1$logFC <- rnorm(length(gr1))
#' gr1$FDR <- rbeta(length(gr1), 1, 8)
#'
#' gr2 <-  GRanges(paste0("chr1:", seq(51, 250, by = 15)))
#' width(gr2) <- 4
#' gr2$logFC <- rnorm(length(gr2))
#' gr2$FDR <- rbeta(length(gr2), 1, 8)
#' grl <- GRangesList(TF1 = gr1, TF2 = gr2)
#' grl <- addDiffStatus(grl)
#'
#' # Using the defaults
#' plotPairwise(grl, var = "logFC")
#'
#' # Density plots on the side panels
#' plotPairwise(
#'   grl, var = "logFC", xside = "density", yside = "density", side_alpha = 0.5
#' )
#'
#' # Turning off side panels, regression line and correlations
#' plotPairwise(
#'   grl, var = "logFC", xside = "none", yside = "none",
#'   rho_geom = "none", line_col = NULL
#' )
#'
#' # Add colours using the status column
#' plotPairwise(grl, var = "logFC", colour = "status") +
#'   scale_fill_manual(values = rep_len(c("blue", "red", "white", "grey"), 8)) +
#'   guides(fill = "none")
#'
#' @importFrom S4Vectors mcols
#' @importFrom rlang '!!' sym
#' @export
plotPairwise <- function(
        x, var, colour = NULL, label = NULL, index = c(1, 2), name_sep = " ",
        xside = c("boxplot", "density", "violin", "none"),
        yside = c("boxplot", "density", "violin", "none"),
        side_panel_width = c(0.3, 0.4), side_alpha = 1,
        xside_axis_pos = "right", yside_axis_label = scales::label_wrap(10),
        line_col = "blue", line_type = 1, line_width = 1,
        rho_geom = c("text", "label", "none"), rho_col = "black",
        rho_size = 4, rho_pos = c(0.05, 0.95),  rho_alpha = 1,
        label_geom = c("label_repel", "label", "text_repel", "text", "none"),
        label_width = 20, label_sep = "; ", label_size = 3.5, label_alpha = 0.7,
        min_d = 1, group_sep = " - ", simplify_equal = TRUE,
        plot_theme = theme_get(), ...
) {
    ## Basic checks
    stopifnot(is(x, "GRangesList"))
    stopifnot(max(index) <= length(x))
    stopifnot(length(x) == 2)
    x <- tryCatch(x[index])
    nm <- names(x)
    stopifnot(length(nm) == 2)
    stopifnot(is(plot_theme, "theme"))

    mc_names <- names(mcols(x[[1]]))
    var <- match.arg(var, mc_names)
    stopifnot(is.numeric(mcols(x[[1]])[[var]]))
    if (!is.null(colour)) colour <- sym(match.arg(colour, mc_names))
    if (!is.null(label)) label <- sym(match.arg(label, mc_names))
    x_lab <- paste(nm[[1]], var, sep = name_sep)
    y_lab <- paste(nm[[2]], var, sep = name_sep)

    xside <- match.arg(xside)
    yside <- match.arg(yside)
    side_panel_width <- rep_len(side_panel_width, 2)
    rho_geom <- match.arg(rho_geom)
    label_geom <- match.arg(label_geom)

    ol <- .makeOLaps(x, var, x_lab, y_lab) ## The basic df for plotting
    if (!is.null(colour)) {
        ol <- .addColourCol(
            x, ol, colour, x_lab, y_lab, group_sep, simplify_equal
        )
    }

    df <- dplyr::filter(ol, !!sym("detected") == "Both Detected")
    p <- ggplot(df, aes(!!sym(x_lab), !!sym(y_lab))) +
        geom_point(aes(colour = {{ colour }}), ...) +
        plot_theme

    if (xside != "none")
        p <- .addXSide(p, ol, x, xside, colour, side_alpha, xside_axis_pos)

    if (yside != "none")
        p <- .addYSide(p, ol, x, yside, colour, side_alpha, yside_axis_label)

    if (xside != "none" | yside != "none") {
        w_x <- side_panel_width[[1]]
        w_y <- side_panel_width[[2]]
        p <- p + theme(
            ggside.panel.scale.x = w_x, ggside.panel.scale.y = w_y,
            axis.title.x = element_text(hjust = 0.5 * (1 - w_y)),
            axis.title.y = element_text(hjust = 0.5 * (1 - w_x))
        )
    }
    if (!is.null(line_col)) {
        p <- p +
            geom_smooth(
                method = "lm", se = FALSE, formula = y ~ x,
                colour = line_col, linetype = line_type, linewidth = line_width
            )
    }
    if (rho_geom != "none")
        p <- .addRho(
            p, ol[[x_lab]], ol[[y_lab]], rho_geom, rho_pos, rho_col, rho_size,
            rho_alpha
        )
    if (!is.null(label) & label_geom != "none") {
        p <- .addLabels(
            x, p, ol, label, colour, x_lab, y_lab, min_d, label_alpha,
            label_size, label_sep, label_geom, label_width
        )
    }
    p + labs(fill = c())

}

#' @importFrom dplyr bind_rows
#' @importFrom S4Vectors mcols
#' @importFrom stringr str_trunc
#' @importFrom ggrepel geom_text_repel geom_label_repel
#' @importFrom rlang '!!'
#' @keywords internal
.addLabels <- function(
        x, p, ol, label, colour, x_lab, y_lab, min_d, .alpha, .size, .sep, .geom,
        .width
) {
    d <- c() # R CMD check
    grp_col <- "detected"
    if (!is.null(colour)) grp_col <- as.character(colour)
    nm <- names(x)
    ol$d <- ol[[x_lab]]^2 + ol[[y_lab]]^2
    split_ol <- split(ol, f = ol[[as.character(grp_col)]])
    split_ol <- lapply(split_ol, dplyr::filter, d == max(d), d > min_d^2)
    split_ol <- split_ol[vapply(split_ol, nrow, integer(1)) == 1]
    ## Now there should be a list of single-row tibbles without NA values
    split_ol <- lapply(
        split_ol,
        function(df) {
            ind <- unlist(df[nm])
            lb <- c(
                unlist(mcols(x[[1]])[[as.character(label)]][ind[[1]]]),
                unlist(mcols(x[[2]])[[as.character(label)]][ind[[2]]])
            )
            lb <- paste(unique(lb), collapse = .sep)
            if (length(lb) == 0) return(NULL)
            df[[as.character(label)]] <- str_trunc(lb, width = .width)
            df
        }
    )
    ol <- bind_rows(split_ol)
    if (nrow(ol) == 0) return(p)
    f <- match.fun(paste0("geom_", .geom))
    p + f(
        aes(label = !!label, colour = {{ colour }}),
        data = ol, alpha = .alpha, size = .size, show.legend = FALSE
    )

}

#' @importFrom stats cor
#' @keywords internal
.addRho <- function(p, xvals, yvals, geom, pos, colour, size, alpha) {
    pos <- rep_len(pos, 2)
    stopifnot(is.numeric(pos))
    rng_x <- range(xvals, na.rm = TRUE)
    rng_y <- range(yvals, na.rm = TRUE)
    cor_x <- min(rng_x) + pos[[1]] * diff(rng_x)
    cor_y <- min(rng_y) + pos[[2]] * diff(rng_y)
    rho <- round(cor(xvals, yvals, use = "pairwise.complete.obs"), 3)
    lab <- paste("rho ==", rho)
    p + annotate(
        geom, x = cor_x, y = cor_y, label = lab, parse = TRUE,
        colour = colour, size = size, alpha = alpha
    )
}

#' @import ggside
#' @importFrom S4Vectors mcols
#' @importFrom rlang '!!' sym
#' @importFrom forcats fct_relabel fct_na_value_to_level
#' @importFrom stringr str_replace_na
#' @keywords internal
.addXSide <- function(p, ol, x, xside, xside_var, alpha, label_side) {
    ## NB: This will be drawn above the plot
    nm <- names(x)
    if (is.null(xside_var)) {
        col <- "detected"
    } else{
        mcol <- as.character(xside_var)
        col <- paste(nm[[2]], mcol, sep = "_")
        vals <- mcols(x[[2]])[[mcol]][ol[[nm[[2]]]]]
        if (is.factor(vals)) {
            vals <- fct_na_value_to_level(vals, "Undetected")
            vals <- fct_relabel(vals, function(x) paste(nm[[2]], x))
        } else {
            vals <- str_replace_na(vals, "Undetected")
            vals <- paste(nm[[2]], vals)
        }
        ol[[col]] <- vals
    }
    stopifnot(col %in% colnames(ol))
    ol <- droplevels(dplyr::filter(ol, !is.na(!!sym(nm[[1]]))))
    if (xside == "density") {
        x_lab <- p$labels$x
        p <- p + geom_xsidedensity(
            aes(x = !!sym(x_lab), y = after_stat(density), fill = !!sym(col)),
            data = ol, colour = NA, alpha = alpha
        )
    }
    if (xside == "boxplot") {
        p <- p + geom_xsideboxplot(
            aes(y = !!sym(col), fill = !!sym(col)), data = ol,
            orientation = "y", alpha = alpha
        ) +
            scale_xsidey_discrete(position = label_side)
    }
    if (xside == "violin") {
        p <- p + geom_xsideviolin(
            aes(y = !!sym(col), fill = !!sym(col)), data = ol,
            orientation = "y", draw_quantiles = 0.5, trim = FALSE,
            alpha = alpha
        ) +
            scale_xsidey_discrete(position = label_side)
    }

    p

}

#' @import ggside
#' @importFrom S4Vectors mcols
#' @importFrom rlang '!!' sym
#' @importFrom forcats fct_relabel fct_na_value_to_level
#' @importFrom stringr str_replace_na
#' @keywords internal
.addYSide <- function(p, ol, x, yside, yside_var, alpha, lab) {
    ## NB: This will be drawn to the right of the plot
    nm <- names(x)
    if (is.null(yside_var)) {
        col <- "detected"
    } else{
        mcol <- as.character(yside_var)
        col <- paste(nm[[1]], mcol, sep = "_")
        vals <- mcols(x[[1]])[[mcol]][ol[[nm[[1]]]]]
        if (is.factor(vals)) {
            vals <- fct_na_value_to_level(vals, "Undetected")
            vals <- fct_relabel(vals, function(x) paste(nm[[1]], x))
        } else {
            vals <- str_replace_na(vals, "Undetected")
            vals <- paste(nm[[1]], vals)
        }
        ol[[col]] <- vals

    }
    stopifnot(col %in% colnames(ol))
    ol <- droplevels(dplyr::filter(ol, !is.na(!!sym(nm[[2]]))))
    if (yside == "density") {
        y_lab <- p$labels$y
        p <- p + geom_ysidedensity(
            aes(y = !!sym(y_lab), x = after_stat(density), fill = !!sym(col)),
            data = ol, orientation = "y", colour = NA, alpha = alpha
        )
    }
    if (yside == "boxplot") {
        p <- p + geom_ysideboxplot(
            aes(x = !!sym(col), fill = !!sym(col)), data = ol,
            orientation = "x", alpha = alpha
        ) +
            scale_ysidex_discrete(labels = lab)
    }
    if (yside == "violin") {
        p <- p + geom_ysideviolin(
            aes(x = !!sym(col), fill = !!sym(col)), data = ol, alpha = alpha,
            orientation = "x", draw_quantiles = 0.5, trim = FALSE
        ) +
            scale_ysidex_discrete(labels = lab)
    }

    p

}


#' @importFrom S4Vectors mcols
#' @importFrom tidyr unite
#' @importFrom rlang '!!' '!!!' syms
#' @importFrom methods is
#' @importFrom forcats fct_relabel fct_cross fct_na_value_to_level
#' @importFrom stringr str_replace_na
#' @importFrom glue glue
#' @keywords internal
.addColourCol <- function(x, ol, colour, x_lab, y_lab, .sep, simplify = TRUE) {
    nm <- names(x)
    col <- as.character(colour)
    ## These columns must be character/factor. They'll be the same type as a GRL
    is_char <- is(mcols(x[[1]])[[col]], "character")
    is_fact <- is(mcols(x[[1]])[[col]], "factor")
    stopifnot(is_char | is_fact)
    ## Now merge as required
    temp <- lapply(nm, function(i) mcols(x[[i]])[[col]][ol[[i]]])
    which_equal <- which(temp[[1]] == temp[[2]])
    pat <- "({nm[[1]]}|{nm[[2]]}).+{.sep}({nm[[1]]}|{nm[[2]]})"
    if (is_char) {
        temp <- lapply(temp, str_replace_na, "Undetected")
        temp <- lapply(seq_along(temp), function(i) paste(nm[[i]], temp[[i]]))
        names(temp) <- nm
        temp_df <- as_tibble(temp)
        temp_df <- unite(
            temp_df, col = !!colour, !!!syms(nm), sep = .sep, remove = FALSE
        )
        if (simplify) {
            temp_df[[col]][which_equal] <- str_replace_all(
                temp_df[[col]][which_equal], glue(pat), "Both"
            )
        }
    } else {
        temp <- lapply(temp, fct_na_value_to_level, "Undetected")
        temp <- lapply(
            seq_along(temp),
            function(i) fct_relabel(temp[[i]], function(x) paste(nm[[i]], x))
        )
        names(temp) <- nm
        temp_df <- as_tibble(temp)
        temp_df[[col]] <- fct_cross(
            temp_df[[nm[[2]]]], temp_df[[nm[[1]]]], sep = .sep
        )
        if (simplify) {
            lv_same <- levels(droplevels(temp_df[[col]][which_equal]))
            lv_same <- str_replace_all(lv_same, glue(pat), "Both")
            levels(temp_df[[col]]) <- c(levels(temp_df[[col]]), lv_same)
            temp_df[[col]][which_equal] <-  str_replace_all(
                temp_df[[col]][which_equal], glue(pat), "Both"
            )
        }
        temp_df <- droplevels(temp_df)
    }
    ol[[col]] <- temp_df[[col]]
    ol
}

#' @importFrom dplyr left_join distinct case_when
#' @importFrom S4Vectors mcols
#' @importFrom rlang '!!!' syms
#' @keywords internal
.makeOLaps <- function(x, var, x_lab, y_lab) {
    nm <- names(x)
    rln <- "many-to-many"
    consensus <- makeConsensus(x)
    ol1 <- as_tibble(findOverlaps(consensus, x[[1]]))
    names(ol1) <- c("queryHits", nm[[1]])
    ol2 <- as_tibble(findOverlaps(consensus, x[[2]]))
    names(ol2) <- c("queryHits", nm[[2]])
    ol <- tibble(queryHits = seq_along(consensus))
    ol <- left_join(
        ol, ol1, by = "queryHits", multiple = "all", relationship = rln
    )
    ol <- left_join(
        ol, ol2, by = "queryHits", multiple = "all", relationship = rln
    )
    ol <- distinct(ol, !!!syms(names(ol)))
    ol$range <- as.character(consensus)[ol$queryHits]

    ol[[x_lab]] <- mcols(x[[1]])[[var]][ol[[nm[[1]]]]]
    ol[[y_lab]] <- mcols(x[[2]])[[var]][ol[[nm[[2]]]]]
    ol$detected <- case_when(
        !is.na(ol[[nm[[1]]]]) & !is.na(ol[[nm[[2]]]]) ~ "Both Detected",
        is.na(ol[[nm[[1]]]]) & !is.na(ol[[nm[[2]]]]) ~ paste(nm[[2]], "Only"),
        !is.na(ol[[nm[[1]]]]) & is.na(ol[[nm[[2]]]]) ~ paste(nm[[1]], "Only"),
        TRUE ~ NA_character_
    )
    lv <- c("Both Detected", paste(nm, "Only"))
    ol$detected <- factor(ol$detected, levels = lv)
    ol
}
