#' @title Perform set operations retaining mcols
#'
#' @description
#' Perform set operations retaining all mcols from the query range
#'
#' @details
#' This extends the methods provided by \link[GenomicRanges]{setdiff},
#' \link[GenomicRanges]{intersect} and \link[GenomicRanges]{union} so that
#' `mcols` from `x` will be returned as part of the output.
#'
#' Where output ranges map back to multiple ranges in `x`, `CompressedList`
#' columns will be returned.
#' By default, these will be simplified if possible, however this behaviour
#' can be disabled by setting `simplify = FALSE`.
#'
#' All columns will be returned which can also be time-consuming.
#' A wise approach is to only provide columns you require as part of the
#' query ranges `x`.
#'
#' If more nuanced approaches are required, the returned columns can be
#' further modified by many functions included in the `plyranges` package,
#' such as `mutate()`.
#'
#' @param x,y GenomicRanges objects
#' @param ignore.strand If set to TRUE, then the strand of x and y is set to
#' "*" prior to any computation.
#' @param simplify logical(1) If TRUE, any List columns will be returned as
#' vectors where possible. This can only occur if single, unique entries are
#' present in all initial elements.
#' @param ... Not used
#'
#' @return
#' A GRanges object with all mcols returned form the original object.
#' If a range obtained by setdiff maps back to two or more ranges in the
#' original set of Ranges, mcols will be returned as
#' \link[IRanges]{CompressedList} columns
#'
#' @examples
#' x <- GRanges("chr1:1-100:+")
#' x$id <- "range1"
#' y <- GRanges(c("chr1:51-60:+", "chr1:21-30:-"))
#' setdiffMC(x, y)
#' setdiffMC(x, y, ignore.strand = TRUE)
#'
#' # The intersection works similarly
#' intersectMC(x, y)
#'
#' # Union may contain ranges not initially in x
#' unionMC(x, y)
#' unionMC(x, y, ignore.strand = TRUE)
#'
#'
#' @importClassesFrom IRanges CompressedList
#' @importFrom GenomicRanges setdiff findOverlaps
#' @importFrom S4Vectors mcols queryHits subjectHits DataFrame
#' @importFrom BiocParallel bplapply bpparam
#'
#' @export
#' @rdname setoptsMC
#' @aliases setdiffMC
setMethod(
  "setdiffMC", c("GRanges", "GRanges"),
  function(x, y, ignore.strand = FALSE, simplify = TRUE, ...) {

    gr <- GenomicRanges::setdiff(x, y, ignore.strand)
    .mapMcols2Ranges(gr, x, ignore.strand, simplify)

  }
)
#' @export
#' @rdname setoptsMC
#' @aliases intersectMC
setMethod(
  "intersectMC", c("GRanges", "GRanges"),
  function(x, y, ignore.strand = FALSE, simplify = TRUE, ...) {

    gr <- GenomicRanges::intersect(x, y, ignore.strand)
    .mapMcols2Ranges(gr, x, ignore.strand, simplify)

  }
)
#' @importClassesFrom GenomicRanges GRangesList
#' @export
#' @rdname setoptsMC
#' @aliases unionMC
setMethod(
  "unionMC", c("GRanges", "GRanges"),
  function(x, y, ignore.strand = FALSE, simplify = TRUE, ...) {

    gr <- GenomicRanges::union(x, y, ignore.strand)
    .mapMcols2Ranges(gr, x, ignore.strand, simplify)

  }
)
#' @export
#' @rdname setoptsMC
#' @aliases setdiffMC
setMethod("setdiffMC", c("ANY", "ANY"), function(x, y, ...) .errNotImp(x, y))
#' @export
#' @rdname setoptsMC
#' @aliases intersectMC
setMethod("intersectMC", c("ANY", "ANY"), function(x, y, ...) .errNotImp(x, y))
#' @export
#' @rdname setoptsMC
#' @aliases unionMC
setMethod("unionMC", c("ANY", "ANY"), function(x, y, ...) .errNotImp(x, y))

#' @importClassesFrom GenomicRanges GRangesList
#' @importFrom IRanges overlapsAny
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors mcols queryHits subjectHits DataFrame
#' @importFrom BiocParallel bplapply bpparam
.mapMcols2Ranges <- function(.gr, .x, .ignore.strand, .simplify) {

  if (ncol(mcols(.x)) == 0) return(.gr)

  ## Treat the ranges differently if there is a match, or no match
  grl <- as.list(split(.gr, f = overlapsAny(.gr, .x)))
  hits <- findOverlaps(grl[["TRUE"]], .x, ignore.strand = .ignore.strand)
  i <- queryHits(hits)
  j <- subjectHits(hits)

  ## If mcols only has one column, a vector will be returned so this handles
  ## the multi-column and single-column case
  DF <- DataFrame(mcols(.x)[j,])
  DF <- setNames(DF, names(mcols(.x)))

  ## Return columns as CompressedList objects if the new ranges map
  ## to multiple ranges in the original object
  if (any(duplicated(i))) {
    DF <- bplapply(
      DF,
      .returnListColumn, i = i, j = j, .simplify = .simplify,
      BPPARAM = bpparam()
    )

  }
  mcols(grl[["TRUE"]]) <- DF

  ## Now replicate the structure of DF exactly for the ranges from y
  n <- length(grl[["FALSE"]])
  if (n > 0) {
    mcols(grl[["FALSE"]]) <- lapply(
      DF,
      function(x) {
        if (is(x, "list_OR_List")) {
          col <- vector("list", n)
        } else {
          col <- rep(NA, n)
        }
        as(col, is(x)[[1]])
      }
    )
  }
  gr <- unlist(GRangesList(grl))
  names(gr) <- c()
  sort(gr)

}

#' @importFrom S4Vectors splitAsList endoapply
.returnListColumn <- function(x, i, j, .simplify) {

  ## x is a vector of any type
  ## i & j are integers denoting vector positions (from queryHits etc)

  out <- splitAsList(x[j], f = as.factor(i))
  names(out) <- c()
  if (is(x, "list_OR_List")) {
    ## Restructure so we have a merged list of the same type as the input
    out <- lapply(out, setNames, nm = c())
    out <- lapply(out, unlist)
    ## Remove explicit NA values as these prevent simplifying later
    out <- lapply(out, na.omit)
    out <- as(out, is(x)[[1]])
  }

  if (.simplify) {
    ## Now make sure only unique entries are returned, and unlist if possible
    ## Also remove NA values
    out <- endoapply(out, unique)
    all_unique <- all(vapply(out, function(x) length(x) == 1, logical(1)))
    if (all_unique) out <- unlist(out)
  }

  out

}

