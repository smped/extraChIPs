#' @title Collapse a GRangesList adding multiple columns from each element
#'
#' @description
#' Make consensus peaks and add individual columns from each original
#' GRangesList element
#'
#' @details
#' Starting with a GRangesList, make a single GRanges object with select columns
#' from each element added to the new object
#'
#' @param x GRangesList
#' @param var columns to map onto the set of consensus peaks
#' @param collapse Any columns specified here will be simplified into a single
#' column
#' @param collapse_sep, String to separate values whe collapsing columns
#' @param name_sep String to separate values when add ing column names
#' @param ... Passed to makeConsensus
#'
#' @return GRanges object
#'
#' @importFrom S4Vectors mcols 'mcols<-'
#' @export
mapGrlCols <- function(
        x, var = NULL, collapse = NULL, collapse_sep = "; ", name_sep = "_", ...
){
    stopifnot(is(x, "GRangesList"))
    consensus <- granges(makeConsensus(x, ...))
    if (is.null(var) & is.null(collapse)) return(consensus)

    mcnames <- .mcolnames(x[[1]])
    if (!is.null(var)) var <- match.arg(var, mcnames, several.ok = TRUE)
    if (!is.null(collapse))
        collapse <- match.arg(collapse, mcnames, several.ok = TRUE)
    var <- unique(c(var, collapse))

    map <- .mapHits(x, consensus)
    map <- .addCols(x, map, var, collapse, collapse_sep, name_sep )

    gr <- consensus[map$consensus_peak]
    df <- dplyr::select(map, -all_of(c("consensus_peak", names(x))))
    mcols(gr) <- as.data.frame(df)

    ## Return any list columns as S4 lists again
    list_cols <- names(which(vapply(df, is, logical(1), "list")))
    mcols(gr)[list_cols] <- lapply(mcols(gr)[list_cols], as, "CompressedList")
    gr

}


#' @importFrom dplyr bind_rows arrange distinct
#' @importFrom tidyr pivot_wider complete
#' @importFrom rlang sym syms '!!' '!!!'
#' @keywords internal
.mapHits <- function(x, consensus) {
    nm <- names(x)
    stopifnot(length(nm) == length(x))
    ol <- lapply(x, findOverlaps, consensus)
    ol <- lapply(ol, as.data.frame)
    ol <- bind_rows(ol, .id = "query")
    ol <- complete(ol, !!sym("subjectHits"))
    ol <- pivot_wider(
        ol, names_from = "query",
        values_from = "queryHits", values_fill = NA, values_fn = list
    )
    ol <- .unnest_by_col(ol, nm, keep_empty = TRUE)
    ol <- arrange(ol, !!sym("subjectHits"))
    ol <- distinct(ol, !!sym("subjectHits"), !!!syms(nm), .keep_all = TRUE)
    dplyr::select(ol, consensus_peak = !!sym("subjectHits"), all_of(nm))
}

#' @importFrom tidyr pivot_longer unnest
#' @importFrom dplyr distinct summarise left_join arrange bind_cols
#' @importFrom rlang ':=' '!!' sym
#' @keywords internal
.addCols <- function(x, map, var, collapse, collapse_sep, name_sep = "_") {
    nm <- names(x)
    if (is.character(collapse_sep)) {
        collapse_sep <- as.list(rep_len(collapse_sep, length(collapse)))
        names(collapse_sep) <- collapse
    }
    if (any(!collapse %in% names(collapse_sep))) stop(
        "All columns being collapsed need a separator provided"
    )
    var_list <- lapply(
        var,
        function(j) { ## j is each column specified in var
            tbl_list <- lapply(
                nm,
                ## i is each elements of x
                function(i) .coerceList(mcols(x[[i]])[[j]])[map[[i]]]
            )
            names(tbl_list) <- paste(nm, j, sep = name_sep)
            tbl <- as_tibble(tbl_list)
            if (j %in% collapse) {
                tbl$consensus_peak <- map$consensus_peak
                tbl <- pivot_longer(
                    tbl, cols = -all_of("consensus_peak"), values_to = j
                )
                tbl <- unnest(tbl, !!sym(j), keep_empty = TRUE)
                tbl <- distinct(tbl, !!sym("consensus_peak"), !!sym(j))
                tbl <- summarise(
                    tbl,
                    "{j}" := paste(
                        sort(!!sym(j)), collapse = collapse_sep[[j]]
                    ),
                    .by = !!sym("consensus_peak")
                )
                tbl <- left_join(
                    dplyr::select(map, !!sym("consensus_peak")), tbl,
                    by = "consensus_peak"
                )
                tbl <- arrange(tbl, !!sym("consensus_peak"))
                tbl <- dplyr::select(tbl, !!sym(j))
            }
            tbl
        }
    )
    df <- bind_cols(var_list)
    bind_cols(map, df)
}


#' @keywords internal
.coerceList <- function(x) {
    ## Not sure shy vctrs::vec_proxy isn't working here
    if (is(x, "List")) x <- as.list(x)
    x
}


#' @importFrom tidyr unnest
.unnest_by_col <- function(data, cols, ...) {
    cols <- match.arg(cols, colnames(data), several.ok = TRUE)
    for (i in cols) data <- unnest(data, all_of(i), ...)
    data
}
