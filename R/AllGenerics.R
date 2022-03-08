setGeneric(
  "colToRanges", function(x, var, ...) standardGeneric("colToRanges")
)

setGeneric(
  "dualFilter", function(x, bg, ref, ...) standardGeneric("dualFilter")
)

setGeneric( "plotAssayRle", function(x, ...) standardGeneric("plotAssayRle"))

setGeneric(
  "setdiffMC",
  function(x, y, ignore.strand = FALSE, simplify = TRUE, ...) {
    standardGeneric("setdiffMC")
  }
)

setGeneric(
  "intersectMC",
  function(x, y, ignore.strand = FALSE, simplify = TRUE, ...) {
    standardGeneric("intersectMC")
  }
)

setGeneric(
  "unionMC",
  function(x, y, ignore.strand = FALSE, simplify = TRUE, ...) {
    standardGeneric("unionMC")
  }
)

setGeneric(
  "reduceMC",
  function(x, ignore.strand = FALSE, simplify = TRUE, ...) {
    standardGeneric("reduceMC")
  }
)

setGeneric(
  "overlapsProp",
  function(x, y, ignore.strand = FALSE, ...) standardGeneric("overlapsProp")
)

setGeneric(
  "partitionRanges",
  function(x, y, ignore.strand = FALSE, suffix = c(".x", ".y"), ...) {
    standardGeneric("partitionRanges")
  }
)

.errNotImp <- function(...){
  args <- list(...)
  cl <- vapply(args, class, character(1))
  msg <- paste("Method not implemented for", paste(cl, collapse = ","))
  message(msg)
}
