#' @export
#' @name bestOverlap
#' @rdname bestOverlap-methods
setGeneric(
    "bestOverlap", function(x, y, ...) standardGeneric("bestOverlap")
)

#' @export
#' @name chopMC
#' @rdname chopMC-methods
setGeneric(
    "chopMC", function(x, simplify = TRUE, ...) standardGeneric("chopMC")
)

#' @export
#' @name colToRanges
#' @rdname colToRanges-methods
setGeneric(
    "colToRanges", function(x, var, ...) standardGeneric("colToRanges")
)

#' @export
#' @name distinctMC
#' @rdname distinctMC-methods
setGeneric(
    "distinctMC",
    function(x, cols, ...) standardGeneric("distinctMC")
)

#' @export
#' @name dualFilter
#' @rdname dualFilter-methods
setGeneric(
    "dualFilter", function(x, bg, ref, ...) standardGeneric("dualFilter")
)

#' @export
#' @name grlToSE
#' @rdname grlToSE-methods
setGeneric("grlToSE", function(x, ...) standardGeneric("grlToSE"))

#' @export
#' @name plotAssayRle
#' @rdname plotAssayRle-methods
setGeneric("plotAssayRle", function(x, ...) standardGeneric("plotAssayRle"))

#' @export
#' @name plotPie
#' @rdname plotPie-methods
setGeneric("plotPie", function(object, ...) standardGeneric("plotPie"))

#' @export
#' @name setoptsMC
#' @rdname setoptsMC-methods
setGeneric(
    "setdiffMC",
    function(x, y, ignore.strand = FALSE, simplify = TRUE, ...) {
        standardGeneric("setdiffMC")
    }
)

#' @export
#' @name setoptsMC
#' @rdname setoptsMC-methods
setGeneric(
    "intersectMC",
    function(x, y, ignore.strand = FALSE, simplify = TRUE, ...) {
        standardGeneric("intersectMC")
    }
)

#' @export
#' @name setoptsMC
#' @rdname setoptsMC-methods
setGeneric(
    "unionMC",
    function(x, y, ignore.strand = FALSE, simplify = TRUE, ...) {
        standardGeneric("unionMC")
    }
)

#' @export
#' @name propOverlap
#' @rdname propOverlap-methods
setGeneric(
    "propOverlap",
    function(x, y, ignore.strand = FALSE, ...) standardGeneric("propOverlap")
)

#' @export
#' @name partitionRanges
#' @rdname partitionRanges-methods
setGeneric(
    "partitionRanges",
    function(
        x, y, y_as_both = TRUE, ignore.strand = FALSE, simplify = TRUE,
        suffix = c(".x", ".y"), ...
    ) {
        standardGeneric("partitionRanges")
    }
)

#' @export
#' @name reduceMC
#' @rdname reduceMC-methods
setGeneric(
    "reduceMC",
    function(x, ignore.strand = FALSE, simplify = TRUE, ...) {
        standardGeneric("reduceMC")
    }
)

.errNotImp <- function(...){
    args <- list(...)
    cl <- vapply(args, class, character(1))
    msg <- paste("Method not implemented for", paste(cl, collapse = ","))
    message(msg)
}
