#' @export
#' @name bestOverlap
#' @rdname bestOverlap-methods
setGeneric("bestOverlap", function(x, y, ...) standardGeneric("bestOverlap"))

#' @export
#' @name colToRanges
#' @rdname colToRanges-methods
setGeneric("colToRanges", function(x, ...) standardGeneric("colToRanges"))

#' @export
#' @name getProfileData
#' @rdname getProfileData-methods
setGeneric(
    "getProfileData", function(x, gr, ...) standardGeneric("getProfileData")
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
#' @name plotSplitDonut
#' @rdname plotSplitDonut-methods
setGeneric(
    "plotSplitDonut", function(object, ...) standardGeneric("plotSplitDonut")
)

#' @export
#' @name setoptsMC
#' @rdname setoptsMC-methods
setGeneric("setdiffMC", function(x, y, ...) standardGeneric("setdiffMC"))

#' @export
#' @name setoptsMC
#' @rdname setoptsMC-methods
setGeneric("intersectMC", function(x, y, ...) standardGeneric("intersectMC"))

#' @export
#' @name setoptsMC
#' @rdname setoptsMC-methods
setGeneric("unionMC", function(x, y, ...) standardGeneric("unionMC"))

#' @export
#' @name propOverlap
#' @rdname propOverlap-methods
setGeneric("propOverlap", function(x, y, ...) standardGeneric("propOverlap"))

#' @export
#' @name partitionRanges
#' @rdname partitionRanges-methods
setGeneric(
    "partitionRanges", function(x, y, ...) standardGeneric("partitionRanges")
)

#' @export
#' @name plotOverlaps
#' @rdname plotOverlaps-methods
setGeneric('plotOverlaps', function(x, ...) standardGeneric('plotOverlaps'))
