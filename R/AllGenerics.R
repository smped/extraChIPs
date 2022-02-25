setGeneric( "plotAssayRle", function(x, ...) standardGeneric("plotAssayRle"))

setGeneric(
  "setdiffMC",
  function(x, y, ignore.strand = FALSE, simplify = TRUE, ...) {
    standardGeneric("setdiffMC")
  }
)

.errNotImp <- function(...){
  args <- list(...)
  cl <- vapply(args, class, character(1))
  msg <- paste("Method not implemented for", paste(cl, collapse = ","))
  message(msg)
}
