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
#' @importFrom grid grid.newpage
f <- function(x, type = c("auto", "venn", "upset"), ...) {

    stopifnot(is(x, "list"))
    stopifnot(length(names(x)) == length(x))
    n <- length(x)
    type <- match.arg(type)
    if (type == "auto") type <- ifelse(n > 3, "upset", "venn")
    x <- lapply(x, unique)

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
        p <- do.call("upset", c(ip, dotArgs))
        return(p)
    }

    if (type == "venn") {
        if (n == 1) p <- .plotSingleVenn(x, ...)
        if (n == 2) p <- .plotDoubleVenn(x, ...)
        if (n == 3) p <- .plotTripleVenn(x, ...)
    }
    p

}
ex <- list(
    x = letters[1:5],
    y = letters[c(6:15, 26)],
    z = letters[c(2, 10:25)]
)
f(ex, type = "upset")

#' @importFrom VennDiagram draw.single.venn
.plotSingleVenn <- function(x, ...) {
    stopifnot(length(x) == 1)
    draw.single.venn(area = length(x[[1]]), category = names(x)[[1]], ...)
}
f(ex[1], fill = "grey")

#' @importFrom VennDiagram draw.pairwise.venn
.plotDoubleVenn <- function(x, ...) {
    stopifnot(length(x) == 2)
    plotArgs <- setNames(lapply(x, length), c("area1", "area2"))
    plotArgs$cross.area <- sum(duplicated(unlist(x)))
    plotArgs$category <- names(x)
    allowed <- c("gList1", "margin", names(formals(draw.pairwise.venn)))
    dotArgs <- list(...)
    dotArgs <- dotArgs[names(dotArgs) %in% allowed]
    do.call("draw.pairwise.venn", c(plotArgs, dotArgs))

}
f(ex[2:3])

#' @importFrom VennDiagram draw.triple.venn
.plotTripleVenn <- function(x, ...) {
    stopifnot(length(x) == 3)
    plotArgs <- setNames(lapply(x, length), paste0("area", 1:3))
    plotArgs$n12 <- sum(duplicated(unlist(x[1:2])))
    plotArgs$n13 <- sum(duplicated(unlist(x[c(1, 3)])))
    plotArgs$n23 <- sum(duplicated(unlist(x[c(2, 3)])))
    plotArgs$n123 <- sum(table(unlist(x)) == 3)
    plotArgs$category <- names(x)
    plotArgs$overrideTriple <- TRUE
    allowed <- c("gList1", "margin", names(formals(draw.triple.venn)))
    dotArgs <- list(...)
    dotArgs <- dotArgs[names(dotArgs) %in% allowed]
    do.call("draw.triple.venn", c(plotArgs, dotArgs))
}
# ex <- list(
#     a = sample(letters, 5),
#     b = sample(letters, 7),
#     c = sample(letters, 10)
# )
f(ex[1:3])


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


