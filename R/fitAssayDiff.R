#' @title Detect Differential ChIP Signal
#'
#' @description Detect differential ChIP signal using one of many approaches
#'
#' @details
#' Starting with a SummarizedExperiment object this function fits either a
#' \link[edgeR]{glmQLFit} model to count data, or the
#' \link[limma:eBayes]{limma-trend} model to logCPM data.
#'
#' If fitting Generalised Linear Models via glmQLFit, options for normalisation
#' are "none", which normalises to library size. Existing library sizes are
#' commonly found in the "totals" column of the colData element and this is
#' attempted by default. All methods provided in \link[edgeR]{calcNormFactors}
#' are also implemented, with the added possibility of normalising within groups
#' instead of across the entire dataset. To enable this, the column with the
#' grouping factor is expected to be in the colData element and is simply
#' called by column name.
#' No normalisation is applied when using the limma-trend model, as this allows
#' for previous normalisation strategies to be performed on the data.
#'
#' Normalising to ChIP Input samples, or using offsets is not yet implemented.
#'
#' Either range-based hypothesis testing is implemented using
#' \link[edgeR]{glmTreat} or \link[limma]{treat}. Setting fc to 1 (or lfc to 0)
#' will default to a point-based null hypothesis, equivalent to either
#' \link[edgeR]{glmQLFTest} (method = "qlf") or
#' \link[limma]{eBayes} (method = "lt").
#'
#' It should also be noted that this is primarily a convenience function and
#' if requiring intermediate output from any setps, then these can be run
#' individually as conventionally specified.
#'
#' @return
#' A SummarizedExperiment object with results set as the `rowData` element.
#' Any existing columns not contained in the differential ChIP results will be
#' retained.
#' Results from testing will contain logCPM, logFC, PValue and the t/F
#' statistic as appropriate, along with an FDR-adjusted p-value
#'
#' @param x a SummarizedExperiment object
#' @param assay The assay to use for analysis
#' @param design The design matrix to use for analysis
#' @param coef The required column from the design matrix
#' @param lib.size The column within the colData element which contains the
#' library size information. If set to NULL, column summaries will be used.
#' @param method the analytic method to be used. Can be 'qlf' which will fit
#' counts using the \link[edgeR]{glmQLFit} strategy , or 'lt' which fits the
#' \link[limma:eBayes]{limma-trend} model on logCPM, or pre-processed logCPM
#' values
#' @param norm The normalisation strategy to use when running the
#' glmQLF models. The value 'none' relies solely on library-size normalisation,
#' and is the default. All methods available in \link[edgeR]{calcNormFactors}
#' are implemented. Ignored when using method = "lt"
#' @param groups character(1) If a column name is supplied here, group-based
#' normalisation will be applied to GLM models treating data in this column
#' as a grouping factor. Ignored when using method = "lt"
#' @param fc,lfc Thresholds passed to \link[limma]{treat} or
#' \link[edgeR]{glmTreat}
#' @param asRanges logical(1). By default, the returned object will be a
#' `SummarizedExperiment` object with the results added to the `rowData`
#' element. Setting `asRanges = TRUE` will only return the GRanges object from
#' this element
#' @param offset If provided will be used as the offset when the DGEList object
#' is created during model fitting
#' @param ... Passed to \link[edgeR]{calcNormFactors}, \link[edgeR]{estimateDisp}
#' and \link[edgeR]{glmQLFit} when method = "qlf".
#' If method = "lt", instead passed to \link[limma]{lmFit}, \link[limma]{treat},
#' \link[limma]{eBayes}
#'
#' @examples
#' nrows <- 200; ncols <- 6
#' counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
#' colnames(counts) <- paste0("Sample_", seq_len(ncols))
#' df <- DataFrame(treat = c("A", "A", "A", "B", "B", "B"))
#' df$treat <- as.factor(df$treat)
#' se <- SummarizedExperiment(
#'   assays = SimpleList(counts = counts), colData = df
#' )
#' X <- model.matrix(~treat, colData(se))
#' se <- fitAssayDiff(se, design = X, lib.size = NULL)
#' rowData(se)
#'
#'
#' @name fitAssayDiff
#' @rdname fitAssayDiff-methods
#' @export
#'
setGeneric(
    "fitAssayDiff", function(x, ...){standardGeneric("fitAssayDiff")}
)
#' @importFrom SummarizedExperiment colData rowData "rowData<-"
#' @importFrom edgeR glmTreat topTags glmQLFTest
#' @importFrom limma eBayes treat topTable topTreat
#' @importFrom stats model.matrix
#' @rdname fitAssayDiff-methods
#' @export
setMethod(
    "fitAssayDiff",
    signature = signature(x = "SummarizedExperiment"),
    function(
        x, assay = "counts", design = NULL, coef = NULL,
        lib.size = "totals", method = c("qlf", "lt"),
        norm = c("none", "TMM", "RLE", "TMMwsp", "upperquartile"),
        groups = NULL, fc = 1, lfc = log2(fc), asRanges = FALSE,
        offset = NULL, ...
    ) {
        method <- match.arg(method)
        norm <- match.arg(norm)
        args <- colnames(colData(x))
        if (is.null(design)) stop("A design matrix must be specified")
        stopifnot(nrow(design) == ncol(x))
        if (is.null(coef)) coef <- colnames(design)[ncol(design)]
        if (is.numeric(coef)) stopifnot(coef <= ncol(design))
        if (is.character(coef)) stopifnot(coef %in% colnames(design))

        if (method == "qlf") {
            ## Only required for GLM fits
            if (!is.null(groups)) groups <- match.arg(groups, args)
            if (!is.null(lib.size)) lib.size <- match.arg(lib.size, args)
            fit <- .se2DGEGLM(
                x, assay, design, lib.size, norm, groups, offset, ...
            )
            if (lfc == 0) {
                fit <- glmQLFTest(fit, coef = coef)
            } else {
                fit <- glmTreat(fit, coef, lfc = lfc)
            }
            res <- topTags(
                fit, n = nrow(x), adjust.method = "none", sort.by = "none"
            )$table
            res <- res[!colnames(res) %in% "unshrunk.logFC"]
        }
        if (method == "lt") {
            fit <- .se2LT(x, assay, design, ...)
            if (lfc == 0) {
                fit <- eBayes(fit, trend = TRUE, ...)
                res <- topTable(
                    fit, coef = coef, number = nrow(x), sort.by = "none",
                    adjust.method = "none"
                )
            } else {
                fit <- treat(fit, lfc = lfc, trend = TRUE, ...)
                res <- topTreat(
                    fit, coef = coef, number = nrow(x), sort.by = "none",
                    adjust.method = "none"
                )
            }
            res <- res[!colnames(res) %in% c("adj.P.Val", "B")]
            colnames(res) <- gsub("AveExpr", "logCPM", colnames(res))
            colnames(res) <- gsub("P.Value", "PValue", colnames(res))
        }
        keep_cols <- setdiff(colnames(rowData(x)), colnames(res))
        orig <- as_tibble(rowData(x))[,keep_cols]
        res[["FDR"]] <- p.adjust(res[["PValue"]], "fdr")
        rowData(x) <- cbind(orig, res)
        if (asRanges) {
            if (is.null(rowRanges(x))) {
                warning(
                    "No ranges found. Results will be returned in the rowData",
                     "element of the original object"
                )
            } else {
                return(rowRanges(x))
            }
        }
        x
    }
)

#' @importFrom SummarizedExperiment assay
#' @importFrom limma lmFit
.se2LT <- function(x, assay, design, ...){
    ## 1. Create an MArrayM/Elist object
    ## 2. Don't normalise
    mat <- assay(x, assay)
    diff <- mat - as.integer(mat)
    if (all(diff == 0)) stop("Expected non-integer values")
    lmFit(mat, design = design, ...)
}

#' @importFrom SummarizedExperiment assay colData
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmQLFit
.se2DGEGLM <- function(x, assay, design, lib.size, norm, groups, offset, ...) {

    ## 1. Create a DGE list
    ## 2. Normalise as requested
    ## 3. Apply the model returning an object of class DGEGLM
    mat <- assay(x, assay)
    if (any(mat < 0)) stop("Counts cannot contain negative values")
    if (anyNA(mat)) stop("Missing values detected")


    col_df <- as_tibble(colData(x), rownames = "colnames")
    ls <- colSums(mat)
    if (!is.null(lib.size)) ls <- col_df[[lib.size]]
    message("Creating DGE list...")
    dge <- DGEList(counts = mat, lib.size = ls, samples = col_df)
    if (!is.null(offset)) dge$offset <- offset
    if (!is.null(groups)) {
        grp_fac <- as.factor(col_df[[groups]])
        dge$samples$group <- grp_fac
        split_df <- split(dge$samples, dge$samples[[groups]])
        list_cols <- lapply(split_df, rownames)
        message("Calculating group-wise normalisation factors...")
        nf <- lapply(
            list_cols,
            function(i) calcNormFactors(dge$counts[,i], method = norm, ...)
        )
        names(nf) <- NULL
        nf <- unlist(nf)
        dge$samples$norm.factors <- nf[colnames(dge)]
    } else {
        message("Calculating experiment-wide normalisation factors...")
        dge <- calcNormFactors(dge, method = norm, ...)
    }
    message("Estimating dispersions...")
    dge <- estimateDisp(dge, design = design, ...)
    message("Running glmQLFit...")
    glmQLFit(dge, design = design, ...)

}

