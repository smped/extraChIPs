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
#' A ggplot2 object able to be customised with colour scales and themes
#'
#' @param df A `data.frame`
#' @param fill The category/column used to fill the slices of the pie charts
#' @param x The second (optional) category/column to place along the x-axis
#' @param y The final (optional) category/column to plce along the y-axis
#' @param width Scale the width of all charts
#' @param show_total logical(1) Show labels on each pie chart with the tally for
#' that complete chart
#' @param label_fill The background colour for tally labels
#' @param label_alpha Transparency for tally labels
#' @param label_size Size of the tally labels. Passed to
#' \link[ggplot2]{geom_label}
#' @param min_p The minimum proportion of the total required for adding labels.
#' Effectively removes labels from pie charts with few members. Alternatively
#' when only one column is specified, categories below this will not be shown
#' around the edge of the plot
#' @param show_category Show category labels around the edge of the plot if only
#' one category/column is specified
#' @param category_size The size of category labels if only one category/column
#' is specified
#' @param category_colour The colour of category labels if only one column is
#' specified
#' @param category_width Width at which category labels will wrap onto a new
#' line
#' @param ... Not used
#'
#' @examples
#' set.seed(200)
#' df <- data.frame(
#'   feature = sample(c("Promoter", "Enhancer", "Intergenic"), 200, replace = TRUE),
#'   TF1 = sample(c("Up", "Down", "Unchanged"), 200, replace = TRUE),
#'   TF2 = sample(c("Up", "Down", "Unchanged"), 200, replace = TRUE)
#' )
#' plotPie(df, fill = "feature")
#' plotPie(df, fill = "feature", x = "TF1")
#' plotPie(df, fill = "feature", x = "TF1", y = "TF2") +
#'  scale_fill_viridis_d() +
#'  theme_bw()
#'
#' @export
plotPie <- function(
  df, fill, x, y, width = 0.8, show_total = TRUE,
  label_fill = "white", label_alpha = 1, label_size = 3, min_p = 0.01,
  show_category = TRUE, category_size = 3, category_colour = "black",
  category_width = 15, ...
) {

  stopifnot(is.data.frame(df))
  if (missing(fill)) stop("The initial category must be defined as 'fill'\n")

  if (missing(x) & missing(y)) {
    p <- .plotSinglePie(
      df = df, fill = fill, width = width, show_total = show_total,
      .lab_fill = label_fill, .lab_alpha = label_alpha, .lab_size = label_size,
      .text_size = category_size, .text_col = category_colour, .min_p = min_p,
      .show_cat = show_category, .text_width = category_width
    )
  }

  if (missing(x) & !missing(y))
    stop("Charts can only be drawn in rows when using two columns.\n")

  if (!missing(x) & missing(y))
    p <- .plotDoublePie(
      df = df, x = x, fill = fill, width = width, show_total = show_total,
      .lab_fill = label_fill, .lab_alpha = label_alpha, .lab_size = label_size,
      .min_p = min_p
    )

  if (!missing(x) & !missing(y))
    p <- .plotTriplePie(
      df = df, x = x, y = y, fill = fill, width = width,
      show_total = show_total, .lab_fill = label_fill, .lab_alpha = label_alpha,
      .lab_size = label_size, .min_p = min_p
    )

  p

}


#' @importFrom dplyr group_by summarise mutate arrange filter
#' @importFrom ggplot2 ggplot aes geom_col geom_text geom_segment geom_label
#' @importFrom ggplot2 coord_polar theme_void theme
#' @importFrom scales comma percent
#' @importFrom rlang '!!' sym
#' @importFrom stringr str_wrap
.plotSinglePie <- function(
  df, fill, width , show_total, .lab_fill, .lab_alpha, .lab_size, .text_size,
  .text_col, .min_p, .show_cat, .text_width
) {

  fill <- fill[[1]]
  stopifnot(fill %in% colnames(df))
  df[[fill]] <- as.factor(df[[fill]])
  y <- lab <- c() # R CMD check error avoidance

  grp_df <- group_by(df, !!sym(fill))
  summ_df <- summarise(grp_df, n = dplyr::n(), .groups = "drop")
  summ_df <- mutate(summ_df, p = n / sum(n))
  summ_df <- arrange(summ_df, desc(!!sym(fill)))
  summ_df$y <- cumsum(summ_df$p)
  summ_df$lab <- paste0(
    summ_df[[fill]], " (",
    percent(summ_df$p, 0.1),
    ")"
  )
  summ_df$lab <- str_wrap(summ_df$lab, width = .text_width)

  p <- ggplot(
    data = summ_df,
    aes(1, p, fill = !!sym(fill))
  ) +
    geom_col(width = width)

  if (.show_cat) {
    p <- p + geom_text(
      aes(x = 1.2 + width / 2, y = y - 0.5*p, label = lab),
      data = dplyr::filter(summ_df, p > .min_p),
      colour = .text_col, size = .text_size
    ) +
    geom_segment(
      aes(
        x = 1 + width / 2, xend = 1.1 + width / 2,
        y = y - 0.5*p, yend = y - 0.5*p
      ),
      data = dplyr::filter(summ_df, p > .min_p)
    )
  }

  if (show_total)
    p <- p + geom_label(
      x = 1 - width / 2, y = 0, label = comma(nrow(df), 1),
      alpha = 0.5, fill = "white", size = .lab_size
    )

  p +
    coord_polar("y", start = 0) +
    theme_void()
}

#' @importFrom dplyr group_by summarise mutate filter ungroup
#' @importFrom tidyr pivot_wider
#' @importFrom tidyselect all_of
#' @importFrom ggplot2 ggplot aes geom_label scale_x_continuous
#' @importFrom ggplot2 coord_equal labs
#' @importFrom scales comma
#' @importFrom rlang '!!' sym
#' @importFrom scatterpie geom_scatterpie
.plotDoublePie <- function(
  df, x, fill, width, show_total, .lab_fill, .lab_alpha, .lab_size, .min_p
) {

  fill <- fill[[1]]
  x <- x[[1]]
  stopifnot(all(c(x, fill) %in% colnames(df) ))
  df[[x]] <- as.factor(df[[x]])
  df[[fill]] <- as.factor(df[[fill]])
  r <- N <- c() # R CMD check error avoidance

  grp_df <- group_by(df, !!sym(x), !!sym(fill))
  summ_df <- summarise(grp_df, n = dplyr::n(), .groups = "drop_last")
  summ_df <- ungroup(mutate(summ_df, N = sum(n)))
  summ_df$r <- summ_df$N / sum(summ_df$N)
  summ_df$r <- 0.5 * summ_df$r / max(summ_df$r) # Set the max as 0.5 always
  wide_df <- pivot_wider(
    summ_df, names_from = all_of(fill), values_from = "n", values_fill = 0
  )
  wide_df$x <- as.integer(wide_df[[x]])

  p <- ggplot() +
    geom_scatterpie(
      aes(x, 1, r = width * r),
      data = wide_df,
      cols = levels(df[[fill]])
    ) +
    coord_equal() +
    scale_x_continuous(
      breaks = seq_along(levels(df[[x]])),
      labels = levels(df[[x]])
    ) +
    labs(x = x, fill = fill) +
    theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  if (show_total)
    p <- p + geom_label(
      aes(x, 1, label = comma(N, 1)),
      data = dplyr::filter(wide_df, N > .min_p * sum(N)),
      size = .lab_size, alpha = .lab_alpha, fill = .lab_fill
    )

  p

}

#' @importFrom dplyr group_by summarise mutate filter ungroup
#' @importFrom tidyr pivot_wider
#' @importFrom tidyselect all_of
#' @importFrom ggplot2 ggplot aes geom_label scale_x_continuous
#' @importFrom ggplot2 coord_equal labs scale_y_continuous
#' @importFrom scales comma
#' @importFrom rlang '!!' sym
#' @importFrom scatterpie geom_scatterpie
.plotTriplePie <- function(
  df, x, y, fill, width, show_total, .lab_fill, .lab_alpha, .lab_size, .min_p
) {

  fill <- fill[[1]]
  x <- x[[1]]
  y <- y[[1]]

  stopifnot(all(c(x, y, fill) %in% colnames(df) ))
  df[[x]] <- as.factor(df[[x]])
  df[[y]] <- as.factor(df[[y]])
  df[[fill]] <- as.factor(df[[fill]])
  r <- N <- c() # R CMD check error avoidance

  grp_df <- group_by(df, !!sym(x), !!sym(y), !!sym(fill))
  summ_df <- summarise(grp_df, n = dplyr::n(), .groups = "drop_last")
  summ_df <- ungroup(mutate(summ_df, N = sum(n)))
  summ_df$r <- summ_df$N / sum(summ_df$N)
  summ_df$r <- 0.5 * summ_df$r / max(summ_df$r) # Set the max as 0.5 always

  wide_df <- pivot_wider(
    summ_df, names_from = all_of(fill), values_from = "n", values_fill = 0
  )
  wide_df$x <- as.integer(wide_df[[x]])
  wide_df$y <- as.integer(wide_df[[y]])

  p <- ggplot() +
    geom_scatterpie(
      aes(x, y, r = width * r),
      data = wide_df,
      cols = levels(df[[fill]])
    ) +
    coord_equal() +
    scale_x_continuous(
      breaks = seq_along(levels(df[[x]])),
      labels = levels(df[[x]])
    ) +
    scale_y_continuous(
      breaks = seq_along(levels(df[[y]])),
      labels = levels(df[[y]])
    ) +
    labs(x = x, y = y, fill = fill)
  if (show_total)
    p <- p + geom_label(
      aes(x, y, label = comma(N, 1)),
      data = dplyr::filter(wide_df, N > .min_p * sum(N)),
      size = .lab_size, alpha = .lab_alpha, fill = .lab_fill
    )

  p

}

