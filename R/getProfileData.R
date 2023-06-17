#' @title Get Profile Data surrounding specified ranges
#'
#' @description Get coverage Profile Data surrounding specified ranges
#'
#' @details
#' This will take all provided ranges and set as identical width ranges,
#' extending by the specified amount both up and downstream of the centre of the
#' provided ranges. By default, the ranges extensions are symmetrical and only
#' the upstream range needs to be specified, however this parameterisation
#' allows for non-symmetrical ranges to be generated.
#'
#' These uniform width ranges will then be used to extract the value contained
#' in the score field from one or more BigWigFiles. Uniform width ranges are
#' then broken into bins of equal width and the average score found within each
#' bin.
#'
#' The binned profiles are returned as a DataFrameList called `profile_data` as
#' a column within the resized GRanges object.
#' Column names in each DataFrame are `score`, `position` and `bp`.
#'
#' If passing a BigWigFileList, profiles will be obtained in series by
#' default. To run in parallel pass a \link[BiocParallel]{MulticoreParam} object
#' to the `BPPARAM` argument.
#'
#' @param x A BigWigFile or BigWiFileList
#' @param gr A GRanges object
#' @param upstream The distance to extend upstream from the centre of each
#' range within `gr`
#' @param downstream The distance to extend downstream from the centre of each
#' range within `gr`
#' @param bins The total number of bins to break the extended ranges into
#' @param mean_mode The method used for calculating the score for each bin.
#' See \link[EnrichedHeatmap]{normalizeToMatrix} for details
#' @param log logical(1) Should the returned values be log2-transformed
#' @param offset Value added to data if log-transforming. Ignored otherwise
#' @param BPPARAM Passed internally to \link[BiocParallel]{bplapply}
#' @param ... Passed to \link[EnrichedHeatmap]{normalizeToMatrix}
#'
#' @return
#' GRanges or GrangesList with column profile_data, as described above
#'
#' @examples
#' bw <- system.file("tests", "test.bw", package = "rtracklayer")
#' gr <- GRanges("chr2:1000")
#' pd <- getProfileData(bw, gr, upstream = 500, bins = 10)
#' pd
#' pd$profile_data
#'
#' @import GenomicRanges
#' @importFrom rtracklayer import.bw BigWigFile BigWigFileList
#' @importFrom EnrichedHeatmap normalizeToMatrix
#' @importFrom tibble as_tibble
#' @importFrom tidyr pivot_longer nest
#' @importFrom tidyselect all_of
#' @importFrom S4Vectors DataFrame mcols
#' @importFrom IRanges SplitDataFrameList
#' @importFrom dplyr left_join
#' @importFrom BiocParallel SerialParam bplapply bpisup bpstart bpstop
#' @importFrom GenomeInfoDb seqlevels keepSeqlevels
#' @rdname getProfileData-methods
#' @export
setMethod(
    "getProfileData",
    signature = signature(x = "BigWigFile", gr = "GenomicRanges"),
    function(
        x, gr, upstream = 2500, downstream = upstream, bins = 100,
        mean_mode = "w0", log = TRUE, offset = 1, ...
    ) {

        stopifnot(upstream > 0 & downstream > 0 & bins > 0)
        ## Manage the seqinfo objects to avoid warnings & errors
        bw_seqlevels <- seqlevels(x)
        gr <- gr[as.character(seqnames(gr)) %in% bw_seqlevels]
        stopifnot(length(gr) > 0)
        bw_seqlevels <- intersect(bw_seqlevels, seqnames(gr))
        gr <- keepSeqlevels(gr, bw_seqlevels)
        ids <- as.character(gr)
        ## Set the bins & resize
        bin_width <- as.integer((upstream + downstream) / bins)
        gr_resize <- resize(gr, width = 1, fix = "center")
        gr_resize <- promoters(gr_resize, upstream, downstream)
        vals <- import.bw(x, which = gr_resize)
        stopifnot("score" %in% colnames(mcols(vals)))
        stopifnot(is.logical(log) & is.numeric(offset))
        if (log) vals$score <- log2(vals$score + offset)

        mat <- normalizeToMatrix(
            signal = vals,
            target = resize(gr_resize, width = 1, fix = "center"),
            extend = (upstream + downstream) / 2, w = bin_width,
            mean_mode = mean_mode, value_column = "score", ...
        )
        mat <- as.matrix(mat)
        rownames(mat) <- ids

        tbl <- as_tibble(mat, rownames = "range")
        tbl <- pivot_longer(
            tbl, cols = all_of(colnames(mat)), names_to = "position",
            values_to = "score"
        )
        tbl[["position"]] <- factor(tbl[["position"]], levels = colnames(mat))
        tbl[["bp"]] <- seq(
            -upstream + bin_width / 2, downstream - bin_width / 2,
            by = bin_width
        )[as.integer(tbl[["position"]])]
        tbl <- nest(tbl, profile_data = all_of(c("score", "position", "bp")))
        gr_tbl <- left_join(as_tibble(gr), tbl, by = "range")
        gr_resize$profile_data <- SplitDataFrameList(
            lapply(gr_tbl$profile_data, DataFrame), compress = TRUE
        )
        gr_resize
    }
)
#' @rdname getProfileData-methods
#' @export
setMethod(
    "getProfileData",
    signature = signature(x = "BigWigFileList", gr = "GenomicRanges"),
    function(
        x, gr, upstream = 2500, downstream = upstream, bins = 100,
        mean_mode = "w0", log = TRUE, offset = 1, BPPARAM = SerialParam(), ...
    ) {
        ## Check BiocParallel is ready to go
        if (!bpisup(BPPARAM)) {
            bpstart(BPPARAM)
            on.exit(bpstop(BPPARAM))
        }
        out <- bplapply(
            x,
            getProfileData,
            gr = gr, upstream = upstream, downstream = downstream, bins = bins,
            mean_mode = mean_mode, log = log, offset = offset, ...,
            BPPARAM = BPPARAM
        )
        as(out, "GRangesList")
    }
)
#' @rdname getProfileData-methods
#' @export
setMethod(
    "getProfileData",
    signature = signature(x = "character", gr = "GenomicRanges"),
    function(
        x, gr, upstream = 2500, downstream = upstream, bins = 100,
        mean_mode = "w0", log = TRUE, offset = 1, ...
    ) {
        stopifnot(all(file.exists(x)))
        if (length(x) == 1) {
            x <- BigWigFile(x)
        } else {
            x <- BigWigFileList(x)
        }
        getProfileData(
            x, gr, upstream, downstream, bins, mean_mode, log, offset, ...
        )
    }
)

