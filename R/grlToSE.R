#' @title Set columns from a GRangesList as Assays in a SummarizedExperiment
#'
#' @description Move one or more columns from a GRangesList elements into
#' assays in a RangesSummarizedEperiment
#'
#' @details
#' Given a GRangesList which would commonly represent multiple samples, reduce
#' any overlapping ranges into a consensus range, setting any metadata columns
#' to be retained as separate assays. These columns may contain values such as
#' coverage, p-values etc.
#'
#' Additional columns can also be placed as rowData columns where the original
#' values are better suited to information about the consensus range rather
#' than the sample (or GRangesList element).
#'
#' Only one value for each range will be retained, and these are chosen using
#' the value provided as the keyvals, taking either the min or max value in this
#' column as the representative range.
#'
#' @param x A GrangesList
#' @param assayCols Columns to move to separate assays
#' @param metaCols Columns to move to mcols within the rowRanges element
#' @param keyvals The value to use when choosing representative values
#' @param by How to choose by keyvals
#' @param ... Passed to \link[GenomicRanges]{reduce}
#' @param ignore.strand logical(1). Whether the strand of the input ranges
#' should be ignored or not.
#'
#' @return
#' A RangedSummarizedExperiment
#'
#' @examples
#' a <- GRanges("chr1:1-10")
#' a$feature <- "Gene"
#' a$p <- 0.1
#' b <- GRanges(c("chr1:6-15", "chr1:15"))
#' b$feature <- c("Gene", "Promoter")
#' b$p <- c(0.5, 0.01)
#' grl <- GRangesList(a = a, b = b)
#' grl
#' se <- grlToSE(
#'   grl, assayCols = "p", metaCols = "feature", keyvals = "p", by = "min"
#' )
#' assay(se, "p")
#' rowRanges(se)
#'
#' @importFrom GenomicRanges reduce
#' @importFrom S4Vectors mcols
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @rdname grlToSE-methods
#' @export
setMethod(
    "grlToSE",
    signature = signature(x = "GRangesList"),
    function(
    x, assayCols = c(), metaCols = c(), keyvals = c(), by = c("min", "max"),
    ..., ignore.strand = FALSE
  ) {

    ## Return an empty object retaining any seqinfo
    if (length(x) == 0) return(.emptySE(x))
    all_mcols <- colnames(mcols(unlist(x)))
    stopifnot(.checkArgsGrlToSe(x, assayCols, metaCols, keyvals, all_mcols))

    by <- match.arg(by)
    assayCols <- intersect(assayCols, all_mcols)
    metaCols <- intersect(metaCols, all_mcols)
    keyvals <- intersect(keyvals, all_mcols)

    nm <- names(x)
    fix_id <- paste0("X", seq_along(x))
    if (is.null(nm)) nm <- fix_id
    if (any(nm == "")) nm[nm == ""] <- fix_id[nm == ""]
    names(x) <- nm

    ## By now we have a named GRL, so get the backbone consensus ranges & sort
    merged_ranges <- reduce(unlist(x), ..., ignore.strand = ignore.strand)
    merged_ranges <- sort(merged_ranges, ignore.strand = ignore.strand)

    assays <- .cols2Assays(
        .merged = merged_ranges, .grl = x, .assayCols = assayCols,
        .keyvals = keyvals, .by = by, .ignore.strand = ignore.strand
    )
    merged_ranges <- .addMcols(merged_ranges, metaCols, x, ignore.strand)

    SummarizedExperiment(assays = assays, rowRanges = merged_ranges)

  }
)

#' @importFrom S4Vectors SimpleList queryHits subjectHits 'mcols<-' mcols
#' @importFrom GenomicRanges findOverlaps GRangesList
#' @importFrom dplyr arrange distinct left_join select across
#' @importFrom rlang '!!!' syms
#' @importFrom tidyr pivot_wider
#' @importFrom tidyselect all_of
.cols2Assays <- function(
    .merged, .grl, .assayCols, .keyvals, .by, .ignore.strand
) {

    ## This needs to be performed for all assayCols together given that values
    ## will be selected based on that in the kay value column
    if (length(.assayCols) == 0) return(SimpleList())
    ## Now move the values to assays using .merged as the backbone
    nm <- names(.grl)
    all_cols <- unique(c(.assayCols, .keyvals))
    ## Return a list with just the values we need
    merged_list <- lapply(
    nm,
    function(i) {
        ## Map each element of the GRList onto the merged backbone in a way that
        ## will return nothing if there's no match
        hits <- findOverlaps(.merged, .grl[[i]], ignore.strand = .ignore.strand)
        gr <- .merged[queryHits(hits)]
        mcols(gr) <- mcols(.grl[[i]][subjectHits(hits)])[all_cols]
        gr$source_in_grl <- i
        gr
    }
    )
    merged_list <- GRangesList(merged_list)
    merged_df <- as_tibble(sort(unlist(merged_list)))

        ## If no key values are given, this will only sort by range
        merged_df <- arrange(merged_df, !!!syms(c("range", .keyvals)))
        # Multiple columns can't be passed to desc so simply reverse here
        if (.by == "max")
        merged_df <- merged_df[seq(nrow(merged_df), 1, by = -1),]
        ## Now by using distinct, we'll effectively have selected the key values
        merged_df <- distinct(
            merged_df,
            across(all_of(c("range", "source_in_grl"))), .keep_all = TRUE
        )

    mats <- lapply(
    .assayCols,
    function(x) {
        df <- pivot_wider(
        merged_df, id_cols = all_of("range"),
        names_from = all_of("source_in_grl"), values_from = all_of(x),
        values_fill = NA
        )
        df <- left_join(tibble(range = as.character(.merged)), df, by = "range")
        df <- dplyr::select(df, all_of(nm))
        stopifnot(length(.merged) == nrow(df)) # Make sure of dims
        stopifnot(length(nm) == ncol(df)) # Make sure of dims
        as.matrix(df)
        }
    )
    names(mats) <- .assayCols
    SimpleList(mats)

}

#' @importFrom rlang '!!' sym
#' @importFrom S4Vectors mcols 'mcols<-' DataFrame
#' @importFrom dplyr group_by summarise across mutate_all
#' @importFrom tidyselect all_of everything
#' @importFrom vctrs vec_proxy
#' @importClassesFrom IRanges CompressedList
.addMcols <- function(.merged, .metaCols, .grl, .ignore.strand) {
    if (length(.metaCols) == 0) return(.merged)
    gr <- unlist(.grl)
    hits <- findOverlaps(.merged, gr, ignore.strand = .ignore.strand)
    df <- data.frame(range = as.character(.merged)[queryHits(hits)])
    df <- cbind(df, as.data.frame(mcols(gr[subjectHits(hits)])[.metaCols]))
    df <- mutate_all(df, vec_proxy) # This handles any AsIs columns
    ## Given that unnesting may behave differently for each required column
    ## unnest within an lapply, then return the vars
    vars <- lapply(
    .metaCols,
    function(x) {
        var_df <- unnest(df[c("range", x)], all_of(x), keep_empty = TRUE)
        var_df <- group_by(var_df, !!sym("range"))
        var <- summarise(var_df, var = list(sort(unique(!!sym(x)))))[["var"]]
        var <- as(var, "CompressedList")
        if (all(vapply(var, length, integer(1)) == 1)) var <- unlist(var)
        var
        }
    )
    names(vars) <- .metaCols
    stopifnot(all(vapply(vars, length, integer(1)) == length(.merged)))
    mcols(.merged) <- DataFrame(vars)
    .merged

}

.emptySE <- function(.x) {
    warning("Returned object will contain no assays or rowData")
    sq <- seqinfo(.x)
    SummarizedExperiment(rowRanges = GRanges(seqinfo = sq))
}

.checkArgsGrlToSe <- function(.x, .assayCols, .metaCols, .keyvals, .all_mcols) {

    ret_val <- TRUE
    msg <- NULL

    if (length(names(.x)) != length(.x))
    msg <- c(msg, "All elements of the x should be given informative names\n")

    if (length(.assayCols) > 0) {
        if (!all(.assayCols %in% .all_mcols)) {
            msg <- c(
                msg,
                paste(
                    "Some columns requested as assays are missing and will be",
                    "ignored.\n"
                )
            )
        }
        if (length(.keyvals) == 0 )
            msg <- c(msg, "Missing keyvals. Values will be chosen at random.\n")

        if (!any(.keyvals %in% .all_mcols)) {
            ret_val <- FALSE
            msg <- c(msg, "No specified 'keyvals'were found.\n")
        }
    }

    if (length(.metaCols) > 0) {
        if (!all(.metaCols %in% .all_mcols)) {
            msg <- c(
                msg,
                paste(
                    "Some columns requested as mcols are missing and will be",
                    "ignored.\n"
                )
            )
        }
    }

    if (any(.metaCols %in% .assayCols)) {
        ret_val <- FALSE
        msg <- c(msg, "Columns for assays and mcols must be distinct.\n")
    }

    if (length(msg) > 0) warning(msg)
    ret_val

}
