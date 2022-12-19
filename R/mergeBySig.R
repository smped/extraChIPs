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
#' @param ... Passed to \link[csaw]{minimalTests}
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
setGeneric(
    "mergeBySig",
    function(x, ...){standardGeneric("mergeBySig")}
)
#' @importFrom S4Vectors DataFrame mcols 'mcols<-'
#' @importFrom GenomicRanges findOverlaps
setMethod(
    "mergeBySig",
    signature = signature(x = "GenomicRanges"),
    function(
        x, df = NULL,
        logfc = "logFC", pval = "P", cpm = "logCPM", inc_cols,
        p_adj_method = "fdr", alpha = 0.05,
        method = c("combine", "best", "minimal"),
        merge_within = 1L, ignore_strand = TRUE,
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
            merged_df <- .ec_combTests(x, ol, df, pval, logfc)
        if (method == "best")
            merged_df <- .ec_bestTest(x, ol, df, pval, logfc, alpha, ret_cols)
        if (method == "minimal")
            merged_df <- .ec_minTest(x, ol, df, pval, logfc, alpha, ...)

        ## Now tidy the standard output from any method
        merged_df <- merged_df[c(n_cols, ret_cols)]
        adj_col <- paste0(pval, "_", p_adj_method)
        merged_df[[adj_col]] <- p.adjust(merged_df[[pval]], p_adj_method)
        mcols(ranges_out) <- merged_df
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


#' @importFrom metapod groupedSimes
#' @importFrom GenomicRanges granges
#' @importFrom S4Vectors queryHits subjectHits DataFrame
#' @importFrom dplyr group_by summarise
#' @importFrom rlang sym '!!'
#' @keywords internal
.ec_combTests <- function(x, ol, df, pval, logfc){
    ## A quick run of microbenchmark to compare csaw::combineTests
    ## showed this is about 25% faster
    ind <- subjectHits(ol)
    grp_df <- as.data.frame(ol)
    simes <- groupedSimes(df[[pval]], ind)
    df$keyval_range <- granges(x)[simes$representative[ind]]
    grp_df <- cbind(grp_df, df[queryHits(ol),])
    grp_df <- group_by(grp_df, subjectHits)
    grp_df$simes_p <- simes$p.value[grp_df$subjectHits]
    merged_df <- summarise(
        grp_df,
        n_windows = dplyr::n(),
        n_up = sum(!!sym(pval) <= simes_p & !!sym(logfc) > 0),
        n_down = sum(!!sym(pval) <=  simes_p & !!sym(logfc) < 0),
        "{pval}" := unique(simes_p),
        .groups = "drop"
    )
    merged_df <- DataFrame(merged_df)
    df <- df[setdiff(names(df), pval)]
    cbind(merged_df, df[simes$representative,])
}

#' @importFrom dplyr group_by mutate summarise across
#' @importFrom tidyselect all_of
#' @importFrom rlang sym '!!'
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomeInfoDb seqinfo
#' @importFrom S4Vectors queryHits DataFrame
#' @keywords internal
.ec_bestTest <- function(x, ol, df, pval, logfc, alpha, ret_cols){
    ## This essentially replicates mergeByCol but using p-values
    ## which is itself csaw::getBestTest(). This is also slightly faster
    grp_df <- as.data.frame(ol)
    ## Needs to be here given the call to 'ret_cols' below
    grp_df[["keyval_range"]] <- as.character(x)[queryHits(ol)]
    grp_df <- cbind(grp_df, df[queryHits(ol),])
    grp_df <- group_by(grp_df, subjectHits)
    grp_df <- mutate(grp_df, holm = p.adjust(!!sym(pval), "holm"))
    merged_df <- summarise(
        grp_df,
        i = which.min(abs(!!sym(pval) - min(!!sym(pval))))[[1]],
        n_windows = dplyr::n(),
        n_up = sum(holm < alpha & !!sym(logfc) > 0),
        n_down = sum(holm < alpha  & !!sym(logfc) < 0),
        across(all_of(ret_cols), function(x) x[i]),
        "{pval}" := holm[i],
        .groups = "drop"
    )
    merged_df <- DataFrame(merged_df)
    sq <- seqinfo(x)
    merged_df$keyval_range <- GRanges(merged_df$keyval_range, seqinfo = sq)
    merged_df
}

#' @importFrom csaw minimalTests
#' @importFrom S4Vectors queryHits DataFrame
#' @importFrom GenomicRanges granges
#' @keywords internal
.ec_minTest <- function(x, ol, df, pval, logfc, alpha, ...) {
    ## Just use the standard csaw function here for simplicity
    ind <- subjectHits(ol)
    min_df <- minimalTests(
        ind, df, pval.col = pval, fc.col = logfc, fc.threshold = alpha, ...
    )
    i <- min_df[["rep.test"]]
    merged_df <- data.frame(
        n_windows = min_df$num.tests,
        n_up = min_df$num.up.logFC,
        n_down = min_df$num.down.logFC
    )
    merged_df <- DataFrame(cbind(merged_df, df[i,]))
    ## Replace the p-values with the holm's one from csaw
    merged_df[[pval]] <- min_df[[pval]]
    merged_df[["keyval_range"]] <- granges(x)[i]
    merged_df
}

