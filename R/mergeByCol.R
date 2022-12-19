#' @title Merge sliding windows using a specified column
#'
#' @description Merge sliding windows using a specified column
#'
#' @details
#' This merges sliding windows using the values in a given column to select
#' representative values for the subsequent merged windows.
#' Values can be chosen from the specified column using any of `min()`,
#' `max()`, `mean()` or `median()`, although `max()` is strongly recommended
#' when specifying values like logCPM.
#' Once a representative range is selected using the specified column, values
#' from columns specified using `inc_cols` are also returned.
#' In addition to these columns, the range from the representative window is
#' returned in the mcols element as a GRanges object in the column
#' `keyval_range`.
#'
#' Merging windows using either the logFC or p-value columns is not implemented.
#'
#' If adjusted p-values are requested an additional column names the same as
#' the initial p-value, but tagged with the adjustment method, will be added.
#' In addition, using the p-value from the selected window, the number of
#' windows with lower p-values are counted by direction and returned in the
#' final object.
#' The selected window will always be counted as up/down regardless of
#' significance as the p-value for this column is taken as the threshold.
#' This is a not dissimilar approach to \link[csaw]{cluster-direction}.
#'
#' If called on a SummarizedExperiment object, the function will be applied to
#' the `rowRanges` element.
#'
#' @param x A GenomicRanges or SummarizedExperiment object
#' @param df A data.frame-like object containing the columns of interest. If
#' not provided, any columns in the mcols() slot will be used.
#' @param col The column to select as representative of the merged ranges
#' @param by The method for selecting representative values
#' @param logfc Column containing logFC values
#' @param pval Column containing p-values
#' @param inc_cols Any additional columns to return. Output will always include
#' columns specified in the arguments `col`, `logfc` and `pval`. Note that
#' values from any additional columns will correspond to the selected range
#' returned in keyval_range
#' @param p_adj_method Any of \link{p.adjust.methods}
#' @param merge_within Merge any ranges within this distance
#' @param ignore_strand Passed internally to \link[GenomicRanges]{reduce} and
#' \link[GenomicRanges]{findOverlaps}
#' @param ... Not used
#'
#' @return
#' A Genomic Ranges object
#'
#' @examples
#' x <- GRanges(c("chr1:1-10", "chr1:6-15", "chr1:51-60"))
#' set.seed(1001)
#' df <- DataFrame(logFC = rnorm(3), logCPM = rnorm(3,8), p = rexp(3, 10))
#' mergeByCol(x, df, col = "logCPM", pval = "p")
#' mcols(x) <- df
#' x
#' mergeByCol(x, col = "logCPM", pval = "p")
#'
#' @name mergeByCol
#' @rdname mergeByCol-methods
#' @export
#'
setGeneric(
    "mergeByCol",
    function(x, ...){standardGeneric("mergeByCol")}
)
#' @importClassesFrom S4Vectors HitsList
#' @importFrom GenomicRanges findOverlaps reduce granges
#' @importFrom GenomeInfoDb seqinfo "seqinfo<-"
#' @importFrom S4Vectors subjectHits queryHits 'mcols<-' mcols
#' @importFrom dplyr group_by summarise n across
#' @importFrom rlang sym '!!'
#' @importFrom stats p.adjust p.adjust.methods
#' @importFrom tidyselect all_of
#' @rdname mergeByCol-methods
#' @export
setMethod(
    "mergeByCol",
    signature = signature(x = "GenomicRanges"),
    function(
        x, df = NULL, col,
        by = c("max", "median", "mean", "min"),
        logfc = "logFC", pval = "P", inc_cols, p_adj_method = "fdr",
        merge_within = 1L, ignore_strand = TRUE, ...
    ) {

        ## Checks & defining the key columns
        if (is.null(df)) df <- mcols(x)
        stopifnot(nrow(df) == length(x))
        df_cols <- colnames(df)
        if(missing(col)) stop("The column containing signal must be specified")
        col <- match.arg(col, df_cols)
        logfc <- match.arg(logfc, df_cols)
        pval <- match.arg(pval, df_cols)
        p_adj_method <- match.arg(p_adj_method, p.adjust.methods)
        f <- match.fun(match.arg(by))
        if (col == logfc | col == pval) stop(
            "Merging not implemented for ", col,
            ". ",
            "Please specify the column containing signal (e.g. logCPM/AveExpr)."
        )

        ## Define the columns to return
        ret_cols <- c("keyval_range", col, logfc, pval)
        ## Can this be rewritten for tidyeval?
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

        ## Merged the data.frame rows
        i <- NULL # Avoid R CMD check issues
        grp_df <- as.data.frame(ol)
        grp_df[["keyval_range"]] <- queryHits(ol)
        grp_df <- cbind(grp_df, df[queryHits(ol),])
        grp_df <- group_by(grp_df, subjectHits)
        merged_df <- summarise(
            grp_df,
            i = which.min(abs(!!sym(col) - f(!!sym(col))))[[1]],
            n_windows = dplyr::n(),
            n_up = sum(!!sym(pval) <= (!!sym(pval))[i] & !!sym(logfc) > 0),
            n_down = sum(!!sym(pval) <= (!!sym(pval))[i] & !!sym(logfc) < 0),
            across(all_of(ret_cols), function(x) x[i]),
            .groups = "drop"
        )

        DF <- DataFrame(merged_df[c(n_cols, ret_cols)])
        adj_col <- paste0(pval, "_", p_adj_method)
        if (p_adj_method != "none")
            DF[[adj_col]] <- p.adjust(DF[[pval]], method = p_adj_method)
        DF[["keyval_range"]] <- granges(x)[DF[["keyval_range"]]]
        seqinfo(DF[["keyval_range"]]) <- seqinfo(x)
        mcols(ranges_out) <- DF
        ranges_out

    }

)
#' @rdname mergeByCol-methods
#' @export
setMethod(
    "mergeByCol",
    signature = signature(x = "RangedSummarizedExperiment"),
    function(
        x, df = NULL, col,
        by = c("max", "median", "mean", "min"),
        logfc = "logFC", pval = "P", inc_cols, p_adj_method = "fdr",
        merge_within = 1L, ignore_strand = FALSE,
        ...
    ) {

        gr <- rowRanges(x)
        mergeByCol(
            gr, df, col, by, logfc, pval, inc_cols, p_adj_method, merge_within,
            ignore_strand, ...
        )

    }

)
