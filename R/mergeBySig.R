#' @title Merge overlapping ranges based on p-values
#'
#' @description Merge overlapping windows using p-values from significance
#' testing
#'
#' @details
#' When using sliding windows to test for differential signal, overlapping
#' windows can be merged based on the significance of results.
#' `mergeBySig()` is a wrapper to the functions \link[csaw]{combineTests},
#' \link[csaw]{getBestTest} and \link[csaw]{minimalTests}, using each
#' function's approach to finding a representative window. The returned object
#' differs from those returned by the original functions in that the
#' description of windows as 'up', 'down' or mixed is omitted and the genomic
#' range corresponding to the representative window is also returned. Column
#' names also correspond to those in the original object.
#'
#' An additional column with adjusted p-values is returned. This column retains
#' the same name as the original but with the suffix '_*' added where the
#' p-value adjustment method is added after the underscore.
#'
#' @param x GenomicRanges object
#' @param df data.frame with results of differential binding analysis performed
#' using a sliding window strategy. If not provided, the columns in the
#' `mcols()` element of `x` will be used
#' @param logfc,pval,cpm Column names for the values holding window specific
#' estimates of change in binding (logfc), overall signal intensity (cpm) and
#' the significance from statistical testing (pval)
#' @param inc_cols (Optional) Character vector of any additional columns in
#' `df` to return
#' @param p_adj_method One of `p.adjust.methods`
#' @param alpha Significance threshold to apply during internal calculations
#' @param method Shorthand versions for which `csaw` strategy to use for
#' merging windows. Choose from 'combine' (\link[csaw]{combineTests}), 'best'
#' (\link[csaw]{getBestTest}) or 'minimal' (\link[csaw]{minimalTests}).
#' @param merge_within Merge any non-overlapping windows within this distance
#' @param ignore_strand Passed internally to \link[GenomicRanges]{reduce} and
#' \link[GenomicRanges]{findOverlaps}
#' @param min_win Only keep merged windows derived from at least this number
#' @param ... Passed to all csaw functions being wrapped
#'
#' @return
#' A GenomicRanges object with overlapping ranges from the original object
#' merged and representative values returned. The range corresponding to the
#' representative values is also returned
#'
#' @examples
#' x <- GRanges(c("chr1:1-10", "chr1:6-15", "chr1:51-60"))
#' set.seed(1001)
#' df <- DataFrame(logFC = rnorm(3), logCPM = rnorm(3,8), p = rexp(3, 10))
#' mcols(x) <- df
#' mergeBySig(x, pval = "p", method = "combine")
#' mergeBySig(x, pval = "p", method = "best")
#' mergeBySig(x, pval = "p", method = "min")
#'
#' @name mergeBySig
#' @rdname mergeBySig-methods
#' @export
setGeneric("mergeBySig", function(x, ...) standardGeneric("mergeBySig"))
#' @importFrom S4Vectors DataFrame mcols mcols<-
#' @import GenomicRanges
#' @rdname mergeBySig-methods
#' @export
setMethod(
    "mergeBySig",
    signature = signature(x = "GenomicRanges"),
    function(
        x, df = NULL,
        logfc = "logFC", pval = "P", cpm = "logCPM", inc_cols,
        p_adj_method = "fdr", alpha = 0.05,
        method = c("combine", "best", "minimal"),
        merge_within = 1L, ignore_strand = TRUE, min_win = 1,
        ...
    ){

        ## Checks & defining the key columns
        if (is.null(df)) df <- mcols(x)
        stopifnot(nrow(df) == length(x))
        df <- DataFrame(df)
        df_cols <- colnames(df)
        logfc <- match.arg(logfc, df_cols)
        pval <- match.arg(pval, df_cols)
        cpm <- match.arg(cpm, df_cols)
        p_adj_method <- match.arg(p_adj_method, p.adjust.methods)
        method <- match.arg(method)

        ## Define the columns to return
        ret_cols <- c("keyval_range", cpm, logfc, pval)
        if (!missing(inc_cols)) {
            inc_cols <- vapply(
                inc_cols, match.arg, character(1), choices = df_cols
            )
            ret_cols <- unique(c(ret_cols, inc_cols))
        }
        ## Always return the columns counting the total, up & down windows
        n_cols <- c("n_windows", "n_up", "n_down")

        ## Merge the ranges and get the map back to the original windows
        ranges_out <- GenomicRanges::reduce(
            x, min.gapwidth = merge_within, ignore.strand = ignore_strand
        )
        ol <- findOverlaps(x, ranges_out, ignore.strand = ignore_strand)

        if (method == "combine")
            merged_df <- .ec_combTests(x, ol, df, pval, logfc, alpha, ...)
        if (method == "best")
            merged_df <- .ec_bestTest(x, ol, df, pval, logfc, alpha, ...)
        if (method == "minimal")
            merged_df <- .ec_minTest(x, ol, df, pval, logfc, alpha, ...)

        ## Now tidy the standard output from any method
        merged_df <- merged_df[c(n_cols, ret_cols)]
        mcols(ranges_out) <- merged_df
        ranges_out <- ranges_out[ranges_out$n_windows >= min_win]
        adj_col <- paste0(pval, "_", p_adj_method)
        vals <- p.adjust(mcols(ranges_out)[[pval]], p_adj_method)
        mcols(ranges_out)[[adj_col]] <- vals

        ranges_out

    }
)
#' @rdname mergeBySig-methods
#' @export
setMethod(
    "mergeBySig",
    signature = signature(x = "RangedSummarizedExperiment"),
    function(
        x, df = NULL,
        logfc = "logFC", pval = "P", cpm = "logCPM", inc_cols,
        p_adj_method = "fdr", alpha = 0.05,
        method = c("combine", "best", "minimal"),
        merge_within = 1L, ignore_strand = TRUE,
        ...
    ) {

        gr <- rowRanges(x)
        mergeBySig(
            gr, df, logfc, pval, cpm, inc_cols, p_adj_method, alpha,
            method, merge_within, ignore_strand, ...
        )

    }

)


#' @import GenomicRanges
#' @importFrom S4Vectors subjectHits DataFrame
#' @importFrom csaw combineTests
#' @keywords internal
.ec_combTests <- function(x, ol, df, pval, logfc, alpha, ...){

    ids <- subjectHits(ol)
    ct <- combineTests(
        ids, df, pval.col = pval, fc.col = logfc, fc.threshold = alpha, ...
    )
    i <- ct[["rep.test"]]
    DF <- DataFrame(
        n_windows = ct[["num.tests"]], n_up = ct[["num.up.logFC"]],
        n_down = ct[["num.down.logFC"]], keyval_range = granges(x)[i]
    )
    DF[[pval]] <- ct[[pval]]
    cols <- setdiff(colnames(df), c(pval))
    rc <- DataFrame(df[i, cols])
    names(rc) <- cols
    cbind(DF, rc)
}

#' @import GenomicRanges
#' @importFrom S4Vectors subjectHits DataFrame
#' @importFrom csaw getBestTest
#' @keywords internal
.ec_bestTest <- function(x, ol, df, pval, logfc, alpha, ...){
    ## This essentially replicates mergeByCol but using p-values
    ## which is itself csaw::getBestTest().
    ids <- subjectHits(ol)
    bt <- getBestTest(
        ids, df, pval.col = pval, fc.col = logfc, fc.threshold = alpha, ...
    )
    i <- bt[["rep.test"]]
    DF <- DataFrame(
        n_windows = bt[["num.tests"]], n_up = bt[["num.up.logFC"]],
        n_down = bt[["num.down.logFC"]], keyval_range = granges(x)[i]
    )
    DF[[pval]] <- bt[[pval]]
    cols <- setdiff(colnames(df), c(pval))
    rc <- DataFrame(df[i, cols])
    names(rc) <- cols
    cbind(DF, rc)
}

#' @importFrom csaw minimalTests
#' @importFrom S4Vectors subjectHits DataFrame
#' @import GenomicRanges
#' @keywords internal
.ec_minTest <- function(x, ol, df, pval, logfc, alpha, ...) {
    ## Just use the standard csaw function here for simplicity
    ids <- subjectHits(ol)
    min_df <- minimalTests(
        ids, df, pval.col = pval, fc.col = logfc, fc.threshold = alpha, ...
    )
    i <- min_df[["rep.test"]]
    merged_df <- DataFrame(
        n_windows = min_df$num.tests, n_up = min_df$num.up.logFC,
        n_down = min_df$num.down.logFC, keyval_range = granges(x)[i]
    )
    merged_df <- cbind(merged_df, df[i,])
    ## Replace the p-values with the holm's one from csaw
    merged_df[[pval]] <- min_df[[pval]]
    merged_df
}

