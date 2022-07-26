#' This function should give the capability to show overlaps for any number of
#' replicates, or a list of items. For n <= 3 **scaled** Venn/Euler diagrams are
#' required, whereas for n > 3 an UpSet plot will be produced.
#' Additionally, an UpSet plot should be optional for n <= 3, so an argument
#' specifying `type` with values c("auto", "venn", "upset") shuold be included.
#' Thus the function requires dispatch for
#' `type = "venn"`
#' - n = 1
#' - n = 2
#' - n = 3
#' `type = "upset"`
#' - n > 1 (given that upset plots are not valid for n = 1)
#'
#' @importFrom methods is
#' @importFrom ComplexUpset upset upset_set_size upset_default_themes
#' @importFrom ggplot2 geom_text aes element_blank
#' @importFrom scales comma
f <- function(x, type = c("auto", "venn", "upset"), ...) {

    stopifnot(is(x, "list"))
    n <- length(x)
    type <- match.arg(type)
    if (type == "auto") type <- ifelse(n > 3, "upset", "venn")

    if (type == "upset") {
        if (n == 1)
            stop("UpSet plots can only be drawn using more than one group")
        ## Setup the df
        all_vals <- as.character(unique(unlist(x)))
        df <- lapply(x, function(i) as.integer(all_vals %in% i))
        df <- as.data.frame(df, row.names = all_vals)
        ip <- list(data = df, intersect = names(df))
        ## Add defaults
        dotArgs <- list(...)
        if (!"set_sizes" %in% names(dotArgs)) {
            dotArgs$set_sizes <- upset_set_size() +
                geom_text(
                    aes(label = comma(..count..)), hjust = 1.15, stat = 'count'
                )
        }
        if (!'themes' %in% dotArgs) {
            dotArgs$themes <- upset_default_themes(panel.grid = element_blank())
        }
        do.call("upset", c(ip, dotArgs))
    }

    if (type == "venn") {
        stop("Not done yet")
        if (n == 1) .plotSingleVenn(x, ...)
        if (n == 2) .plotDoubleVenn(x, ...)
        if (n == 3) .plotTripleVenn(x, ...)
    }

}
f(ex, type = "upset")

#' @importFrom methods is
#' @importFrom tibble as_tibble
#' @importFrom dplyr bind_rows summarise group_by arrange
#' @importFrom stringr str_count
#' @importFrom forcats fct_inorder
.calculateAreas <- function(data, .sep = "-") {

    ## Limit this to 3way interactions. Anything bigger is better as an UpSet plot
    stopifnot(is(data, "list"))
    stopifnot(all(vapply(data, is.character, logical(1))))

    tbl <- bind_rows(lapply(data, as_tibble), .id = "group")
    grp_tbl <- group_by(tbl, value)
    summ_tbl <- summarise(
        grp_tbl,
        group = paste(unique(sort(group)), collapse = .sep),
        .groups = 'drop'
    )
    summ_tbl$n <- str_count(summ_tbl$group, .sep)
    summ_tbl <- arrange(summ_tbl, n, group)
    summ_tbl$group <- fct_inorder(summ_tbl$group)

    out <- split(summ_tbl, summ_tbl$group)
    lapply(out, nrow)

}
ex <- list(
    x = letters[1:10],
    y = letters[c(6:15, 26)],
    z = letters[c(2, 10:25)]
)
.calculateAreas(ex)


