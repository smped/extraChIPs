#' @title Plot Overlaps Between List Elements
#'
#' @description Plot Overlaps between list elements as an upset or Venn diagram
#'
#' @details
#' This function should give the capability to show overlaps for any number of
#' replicates or groups, or a list of items such as gene names.
#' For n = 2, a scaled Venn Diagram will be produced, however no scaling is
#' implemented for n = 3
#'
#' UpSet plots are possible for any lists with length > 1, and are the only
#' implemented possibility for lists > 3.
#'
#' If the input is a `GRangesList` an additional boxplot can be requested
#' using any numeric column within the existing `mcols()` element.
#' Values will be summarised across all elements using the requested function
#' and the boxplot will be included as an upper panel above the intersections
#'
#' @return
#' Either a VennDiagram (i.e. grid) object, or a ComplexUpset plot
#'
#' @param x GRangesList of S3 list able to be coerced to character vectors
#' @param type The type of plot to be produced
#' @param var Column to summarised as a boxplot in an upper panel
#' (UpSet plot only)
#' @param f Summarisation function. Must return a single value from any
#' numeric vector
#' @param set_col Colours to be assigned to each set
#' @param ... Passed to \link[VennDiagram]{draw.pairwise.venn} (or
#' `draw.single/triple.venn`) for Venn Diagrams, and to
#' \link[SimpleUpset]{simpleUpSet} for UpSet plots
#' @param ignore.strand Passed to \link[GenomicRanges]{reduce}
#' @param merge_within Passed to \link{makeConsensus}
#' @param label_size Text size for set and intersection labels. Passed
#' internally to `geom_text(size = label_size)`
#' @param hj_sets Horizontal adjustment of set size labels
#' @param vj_intersect Vertical adjustment of intersection size labels
#' @param exp_sets X-axis expansion for set size panel
#' @param exp_intersect Y-axis expansion for intersections size panel
#'
#' @examples
#' ## Examples using a list of character vectors
#' ex <- list(
#'   x = letters[1:5], y = letters[c(6:15, 26)], z = letters[c(2, 10:25)]
#' )
#' plotOverlaps(ex, type = "upset")
#' plotOverlaps(ex, type = "venn", set_col = 1:3, alpha = 0.3)
#' plotOverlaps(ex[1:2])
#'
#' ## GRangesList object will produce a boxplot of summarised values in the
#' ## upper panel
#' data("peaks")
#' grl <- peaks[1:3]
#' names(grl) <- gsub("_peaks.+", "", names(grl))
#' plotOverlaps(grl, type = 'upset', var = 'score', f = 'max')
#'
#' ## If only two samples are present, a VennDiagram will be produced
#' plotOverlaps(grl[1:2], set_col = c("green", "blue"), cex = 1.5, cat.cex = 1.5)
#'
#' @import GenomicRanges
#' @importFrom S4Vectors endoapply mcols
#' @importFrom dplyr bind_cols
#' @importFrom rlang list2 := !! sym
#' @import ggplot2
#' @rdname plotOverlaps-methods
#' @aliases plotOverlaps
#' @export
setMethod(
  "plotOverlaps",
  signature = "GRangesList",
  function(
    x, type = c("auto", "venn", "upset"), var = NULL,
    f = c("mean", "median", "max", "min", "sd"),
    merge_within = 1L, ignore.strand = TRUE, set_col = NULL, ...,
    label_size = 3.5,
    hj_sets = 1.15, exp_sets = 0.2,
    vj_intersect = - 0.5, exp_intersect = 0.1
  ) {

    stopifnot(methods::is(x, "GRangesList"))
    nm <- names(x)
    n <- length(x)
    stopifnot(length(nm) == n)
    type <- match.arg(type)
    if (type == "auto") type <- ifelse(n > 3, "upset", "venn")

    if (!is.null(var)) var <- match.arg(var[[1]], .mcolnames(x[[1]]))

    # Collapse as required
    gr <- makeConsensus(
      x, var = var, merge_within = merge_within, ignore.strand = ignore.strand
    )

    if (is.null(var) | type == "venn") {

      ## Form a character list & plot
      l <- lapply(nm, function(x) as.character(gr)[mcols(gr)[[x]]])
      names(l) <- nm
      plotOverlaps(
        l, type, set_col = set_col, sz_sets = label_size,
        hj_sets = hj_sets, exp_sets = exp_sets,
        exp_intersect = exp_intersect, vj_intersect = vj_intersect, ...
      )

    } else {

      if (n == 1) stop("UpSet plots can only be drawn using more than one group")

      if (!is.numeric(mcols(x[[1]])[[var]]))
        stop(var, " must contain numeric values")

      ## Setup the df
      tbl <- as_tibble(gr)
      f <- match.arg(f)
      f <- match.fun(f)
      if (methods::is(tbl[[var]], "list"))
        tbl[[var]] <- vapply(tbl[[var]], f, numeric(1))

      p <- .makeUpSet(
        tbl, nm, var, set_col, label_size, hj_sets, exp_sets, exp_intersect,
        vj_intersect, ...
      )
      return(p)

    }

  }
)

#'
#' @rdname plotOverlaps-methods
#' @aliases plotOverlaps
#' @export
setMethod(
  "plotOverlaps",
  signature = "list",
  function(
    x, type = c("auto", "venn", "upset"), set_col = NULL, ..., label_size = 3.5,
    hj_sets = 1.15, exp_sets = 0.2, vj_intersect = - 0.5, exp_intersect = 0.1
  ) {

    stopifnot(length(names(x)) == length(x))
    n <- length(x)
    nm <- names(x)
    type <- match.arg(type)
    if (type == "auto") type <- ifelse(n > 3, "upset", "venn")
    x <- lapply(x, as.character)
    x <- lapply(x, unique)

    if (type == "upset") {

      if (n == 1)
        stop("UpSet plots can only be drawn using more than one group")

      ## Setup the df
      all_vals <- unique(unlist(x))
      df <- lapply(x, function(i) all_vals %in% i)

      ## Ensure colnames are respected
      df <- as.data.frame(df)
      names(df) <- nm

      p <- .makeUpSet(
        df, sets = nm, var = NULL, set_col, label_size, hj_sets, exp_sets,
        exp_intersect, vj_intersect, ...
      )
      return(p)
    }

    if (type == "venn") {
      grid::grid.newpage()
      if (n == 1) p <- .plotSingleVenn(x, fill = set_col, ...)
      if (n == 2) p <- .plotDoubleVenn(x, fill = set_col, ...)
      if (n == 3) p <- .plotTripleVenn(x, fill = set_col, ...)
    }
    invisible(p)

  }
)

#' @keywords internal
.makeUpSet <- function(
    x, sets, var, set_col, label_size, hj_sets, exp_sets,
    exp_intersect, vj_intersect, ...
) {

  if (!requireNamespace('SimpleUpset', quietly = TRUE))
    stop("Please install 'SimpleUpset' to use this function.")

  ## Now use SimpleUpset
  args <- list(x = x, sets = sets)
  if (is.null(set_col)) {
    args$set_layers <- SimpleUpset::default_set_layers(
      hjust = hj_sets, label_size = label_size, expand = c(exp_sets, 0)
    )
  } else {
    args$set_layers <- SimpleUpset::default_set_layers(
      fill = "set", hjust = hj_sets, label_size = label_size,
      expand = c(exp_sets, 0),
      scale_fill_manual(values = rep_len(set_col, length(sets))),
      guides(fill = guide_none())
    )
  }
  args$intersect_layers = SimpleUpset::default_intersect_layers(
    expand = c(0, exp_intersect), label_size = label_size,
    vjust = vj_intersect
  )
  if (!is.null(var)) {
    var <- match.arg(var, colnames(x))
    args$annotation <- list(aes(y = !!sym(var)), geom_boxplot())
  }
  args <- c(args, list(...))
  do.call(SimpleUpset::simpleUpSet, args)

}

.plotSingleVenn <- function(x, ...) {
  if (!requireNamespace('VennDiagram', quietly = TRUE))
    stop("Please install 'VennDiagram' to use this function.")
  stopifnot(length(x) == 1)
  VennDiagram::draw.single.venn(
      area = length(x[[1]]), category = names(x)[[1]], ...
  )
}

.plotDoubleVenn <- function(x, ...) {
  if (!requireNamespace('VennDiagram', quietly = TRUE))
    stop("Please install 'VennDiagram' to use this function.")
  stopifnot(length(x) == 2)
  plotArgs <- setNames(lapply(x, length), c("area1", "area2"))
  plotArgs$cross.area <- sum(duplicated(unlist(x)))
  plotArgs$category <- names(x)
  vd_formals <- names(formals(VennDiagram::draw.pairwise.venn))
  allowed <- c("gList1", "margin", vd_formals)
  dotArgs <- list(...)
  dotArgs <- dotArgs[names(dotArgs) %in% allowed]
  do.call(VennDiagram::draw.pairwise.venn, c(plotArgs, dotArgs))

}

.plotTripleVenn <- function(x, ...) {
  if (!requireNamespace('VennDiagram', quietly = TRUE))
    stop("Please install 'VennDiagram' to use this function.")
  stopifnot(length(x) == 3)
  plotArgs <- setNames(lapply(x, length), paste0("area", seq_len(3)))
  plotArgs$n12 <- sum(duplicated(unlist(x[c(1, 2)])))
  plotArgs$n13 <- sum(duplicated(unlist(x[c(1, 3)])))
  plotArgs$n23 <- sum(duplicated(unlist(x[c(2, 3)])))
  plotArgs$n123 <- sum(table(unlist(x)) == 3)
  plotArgs$category <- names(x)
  plotArgs$overrideTriple <- TRUE
  vd_formals <- names(formals(VennDiagram::draw.triple.venn))
  allowed <- c("gList1", "margin", vd_formals)
  dotArgs <- list(...)
  dotArgs <- dotArgs[names(dotArgs) %in% allowed]
  do.call(VennDiagram::draw.triple.venn, c(plotArgs, dotArgs))
}


#' @importFrom S4Vectors mcols
#' @keywords internal
.mcolnames <- function(x) {colnames(mcols(x))}
