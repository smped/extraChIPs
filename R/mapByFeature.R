#' @title Map Genomic Ranges to genes using defined features
#'
#' @description Map Genomic Ranges to genes using defined regulatory features
#'
#' @details
#' This function is able to utilise feature-level information and long-range
#' interactions to enable better mapping of regions to genes. If provided, this
#' essentially maps from ranges to genes using the regulatory features as a
#' framework. The following sequential strategy is used:
#'
#' \enumerate{
#'   \item Ranges overlapping a promoter are assigned to that gene
#'   \item Ranges overlapping an enhancer are assigned to **all genes** within
#'   a specified distance
#'   \item Ranges overlapping a long-range interaction are assigned to all genes
#'   connected by the interaction
#'   \item Ranges with no gene assignment from the previous steps are assigned
#'   to *all overlapping genes* or the nearest gene within a specified distance
#' }
#'
#' If information is missing for one of these steps, the algorithm will simply
#' proceed to the next step. If no promoter, enhancer or interaction data is
#' provided, all ranges will be simply mapped by step 4.
#' Ranges can be mapped by any or all of the first three steps, but step 4 is
#' mutually exclusive with the first 3 steps.
#'
#' Distances between each set of features and the query range can be
#' individually specified by modifying the `gr2prom`, `gr2enh`, `gr2gi` or
#' `gr2gene` parameters. Distances between features and genes can also be set
#' using the parameters `prom2gene`, `enh2gene` and `gi2gene`.
#'
#' Additionally, if previously defined mappings are included with any of the
#' `prom`, `enh` or `gi` objects, this will be used in preference to any
#' obtained from the `genes` object.
#'
#' @return
#' A GRanges object with added mcols as specified
#'
#' @param gr GRanges object with query ranges to be mapped to genes
#' @param genes GRanges object containing genes (or any other nominal feature)
#' to be assigned
#' @param prom GRanges object defining promoters
#' @param enh GRanges object defining Enhancers
#' @param gi GInteractions object defining interactions. Mappings from
#' interactions to genes should be performed as a separate prior step.
#' @param cols Column names to be assigned as mcols in the output. Columns
#' must be minimally present in `genes`. If all requested columns are found in
#' any of prom, enh or gi, these pre-existing mappings will be preferentially
#' used. Any columns not found in utilised reference objects will be ignored.
#' @param gr2prom The maximum permissible distance between a query range and
#' any ranges defined as promoters
#' @param gr2enh The maximum permissible distance between a query range and
#' any ranges defined as enhancers
#' @param gr2gi The maximum permissible distance between a query range and
#' any ranges defined as GInteraction anchors
#' @param prom2gene The maximum permissible distance between a range provided
#' in `prom` and a gene
#' @param enh2gene The maximum permissible distance between a range provided
#' in `enh` and a gene
#' @param gi2gene The maximum permissible distance between a GInteractions
#' anchor (provided in `gi`) and a gene
#' @param gr2gene The maximum permissible distance between a query range and
#' genes (for ranges not otherwise mapped)
#' @param ... Passed to findOverlaps and overlapsAny internally
#'
#' @examples
#' ## Define some genes
#' genes <- GRanges(c("chr1:2-10:*", "chr1:25-30:-", "chr1:31-40:+"))
#' genes$gene_id <- paste0("gene", seq_along(genes))
#' genes
#' ## Add a promoter for each gene
#' prom <- promoters(genes, upstream = 1, downstream = 1)
#' prom
#' ## Some ranges to map
#' gr <- GRanges(paste0("chr1:", seq(0, 60, by = 15)))
#' gr
#'
#' ## Map so that any gene within 25bp of the range is assigned
#' mapByFeature(gr, genes, gr2gene = 25)
#'
#' ## Now use promoters to be more accurate in the gene assignment
#' ## Given that the first range overlaps the promoter of gene1, this is a
#' ## more targetted approach. Similarly for the third range
#' mapByFeature(gr, genes, prom, gr2gene = 25)
#'
#' @importFrom S4Vectors mcols subjectHits mcols<-
#' @importClassesFrom IRanges CompressedList
#' @importFrom IRanges overlapsAny
#' @importFrom dplyr bind_rows distinct across left_join
#' @importFrom tidyr chop
#' @importFrom tidyselect all_of
#' @importFrom methods as
#' @importFrom stats setNames
#' @import GenomicRanges
#'
#' @export
mapByFeature <- function(
        gr, genes, prom, enh, gi, cols = c("gene_id", "gene_name", "symbol"),
        gr2prom = 0, gr2enh = 0, gr2gi = 0, gr2gene = 1e5,
        prom2gene = 0, enh2gene = 1e5, gi2gene = 0, ...
) {

    ## Add some checks here
    checkArgs <- .checkMappingArgs(
        .gr = gr, .genes = genes, .prom = prom, .enh = enh, .gi = gi,
        .cols = cols, .distArgs = list(
            gr2prom, gr2enh, gr2gi, gr2gene, prom2gene, enh2gene, gi2gene
        )
    )
    if (!checkArgs) stop("Incorrectly specified arguments")

    ## Perform mapping in a sequential manner as described in the algorithm
    ## First identify the mappings from gr to features
    ol_prom <- ol_enh <- ol_gi <- rep(FALSE, length(gr))
    if (!missing(prom)) ol_prom <- overlapsAny(gr, prom, maxgap = gr2prom, ...)
    if (!missing(enh)) ol_enh <- overlapsAny(gr, enh, maxgap = gr2enh, ...)
    if (!missing(gi)) ol_gi <- overlapsAny(gr, gi, maxgap = gr2gi, ...)
    ol_none <- !ol_prom & !ol_enh & !ol_gi

    ## 1. Map promoters
    prom_tbl <- .mapFeatures(
        gr[ol_prom,], prom, genes, cols, gr2prom, prom2gene, ...
    )
    ## 2. Map enhancers
    enh_tbl <- .mapFeatures(
        gr[ol_enh], enh, genes, cols, gr2enh, enh2gene, ...
    )
    ## 3. Add gi mappings
    gi_tbl <- .mapGi(gr[ol_gi], gi, genes, cols,  gr2gi, gi2gene, ...)
    ## 4. Map any unmapped regions
    misc_tbl <- .mapWithin(gr[ol_none], genes, cols, gr2gene, ...)

    ## Return the complete set of mappings & keep them as unique mappings
    mapped_tbl <- bind_rows(
        list(prom_tbl, enh_tbl, bind_rows(gi_tbl), misc_tbl)
    )
    cols <- intersect(cols, names(mapped_tbl))
    mapped_tbl <- distinct(
        mapped_tbl, across(all_of(c("range", cols))), .keep_all = TRUE
    )
    mapped_tbl <- chop(mapped_tbl, all_of(cols))

    ## Now map back onto the original object using tidy methods
    gr_tbl <- tibble(range = as.character(gr))
    gr_tbl <- left_join(gr_tbl, mapped_tbl, by = "range")
    gr_list <- as.list(gr_tbl)
    list_cols <- vapply(gr_tbl, is.list, logical(1))
    gr_list[list_cols] <- lapply(
        ## The new vctrs-based types are problematic. Coerce via as.list()
        gr_list[list_cols], function(x) as(as.list(x), "CompressedList")
    )

    out <- granges(gr)
    mcols(out) <-  cbind(mcols(gr)[!.mcolnames(gr) %in% cols], gr_list[cols])
    out

}

#' @title Map ranges to genes using features as an anchor
#' @param .gr The ranges to map onto
#' @param .feat Features to use for mapping
#' @param .genes GRanges object containing gene-level information
#' @param .cols The columns from `.genes` to map onto `.gr`
#' @param .gr2feat The maximum distance between ranges and features
#' @param .feat2gene The maximum distance between features & genes
#' @param ... Passed to findOverlaps and subsetByOverlaps
#' @return A data.frame
#' @importFrom S4Vectors mcols
#' @importFrom dplyr inner_join mutate_all
#' @importFrom vctrs vec_proxy
#' @importFrom tidyr unnest
#' @importFrom tidyselect everything
#' @import GenomicRanges
#' @importFrom IRanges subsetByOverlaps
#' @keywords internal
.mapFeatures <- function(.gr, .feat, .genes, .cols, .gr2feat, .feat2gene, ...) {

    ## Restrict the data to just the overlapping subsets
    if (missing(.feat)) return(NULL)
    if (is.null(.feat)) return(NULL)
    .gr <- subsetByOverlaps(.gr, .feat, maxgap = .gr2feat)
    if (length(.gr) == 0) return(NULL)
    .feat <- subsetByOverlaps(.feat, .gr, maxgap = .gr2feat, ...)
    if (length(.feat) == 0) return(NULL)

    ## Given that we may pull the data directly from features, check both
    if (!missing(.genes)) {
        .cols <- intersect(.cols, names(mcols(.genes)))
    } else {
        .cols <- intersect(.cols, names(mcols(.feat)))
    }
    if (length(.cols) == 0) stop("Could not find requested column")

    if (!all(.cols %in% names(mcols(.feat)))) {
        ## Use the genes object here. Otherwise it is ignored
        ## Map to genes if any requested column is absent from existing mappings
        .genes <- subsetByOverlaps(.genes, .feat, maxgap = .feat2gene, ...)
        if (length(.genes) == 0) return(NULL) # No point proceeding from here
        feat_df <- findOverlaps(.feat, .genes, maxgap = .feat2gene, ...)
        feat_df <- as.data.frame(feat_df)
        .genes$subjectHits <- seq_along(.genes)
        gene_df <- as.data.frame(mcols(.genes)[c("subjectHits", .cols)])
        feat_df <- inner_join(feat_df, gene_df, by = "subjectHits")
        feat_df <- feat_df[c("queryHits", .cols)]
        ## Rename for mapping onto the ranges
        names(feat_df) <- gsub("queryHits", "subjectHits", names(feat_df))
    } else {
        ## Setup as the same format otherwise
        feat_df <- mcols(.feat)[.cols]
        feat_df[["subjectHits"]] <- seq_along(.feat)
        feat_df <- as.data.frame(feat_df)
    }

    ## Deal with any list columns
    feat_df <- mutate_all(feat_df, vec_proxy) # Remove any AsIs attributes!!!
    feat_df <- unnest(feat_df, everything())

    ## Map ranges to features
    range_to_feat <- findOverlaps(.gr, .feat, maxgap = .gr2feat, ...)
    range_to_feat <- as.data.frame(range_to_feat)
    range_df <- data.frame(
        range = as.character(.gr), queryHits = seq_along(.gr)
    )
    ## This will drop any additional columns. They can be replaced back in the
    ## parent function
    range_df <- inner_join(
        range_df, range_to_feat, by = "queryHits", multiple = "all"
    )
    range_df <- inner_join(
        range_df, feat_df, by = "subjectHits", multiple = "all"
    )
    range_df <- range_df[c("range", .cols)]

    ## Return the tibble
    range_df

}

#' @title Map ranges to genes via Interactions
#' @param .gr The ranges to map onto
#' @param .gi GInteractions object
#' @param .genes GRanges object containing gene-level information
#' @param .cols The columns from `.genes` to map onto `.gr`
#' @param .gr2gi The maximum distance between ranges and anchors
#' @param .gi2gene The maximum distance between anchors & genes
#' @param ... Passed to findOverlaps
#' @return data.frame of mapped ranges
#' @importFrom S4Vectors mcols
#' @importFrom vctrs vec_proxy
#' @importFrom InteractionSet anchors
#' @importFrom IRanges subsetByOverlaps
#' @importFrom dplyr inner_join mutate_all
#' @importFrom tidyr unnest
#' @importFrom tidyselect everything
#' @import GenomicRanges
#' @keywords internal
.mapGi <- function(.gr, .gi, .genes, .cols, .gr2gi, .gi2gene , ...) {

    ## Restrict the data to just the overlapping subsets
    if (missing(.gi)) return(NULL)
    if (is.null(.gi)) return(NULL)
    .gr <- subsetByOverlaps(.gr, .gi, maxgap = .gr2gi)
    if (length(.gr) == 0) return(NULL)
    .gi <- subsetByOverlaps(.gi, .gr, maxgap = .gr2gi, ...)
    if (length(.gi) == 0) return(NULL)

    if (!missing(.genes)) {
        .cols <- intersect(.cols, names(mcols(.genes)))
    } else {
        .cols <- intersect(.cols, names(mcols(.gi)))
    }
    if (length(.cols) == 0) stop("Could not find requested column")

    if (!all(.cols %in% names(mcols(.gi)))) {
        ## If any requested column is not in the interactions object, use the
        ## genes. Just add these as columns then run as per normal
        .genes <- subsetByOverlaps(.genes, .gi, maxgap = .gi2gene)
        if (length(.genes) == 0) return(NULL) # Exit if no genes can be mapped
        anchors_to_gene <- lapply(
            anchors(.gi),
            .mapWithin,
            .genes = .genes, .cols = .cols, .within = .gi2gene
        ) # This returns the anchor range, not .gr
        anchors_to_gene <- bind_rows(anchors_to_gene)
        ## Handle any list columns in case they exist
        anchors_to_gene <- mutate_all(anchors_to_gene, vec_proxy)
        anchors_to_gene <- unnest(anchors_to_gene, everything())
        ## Map back to the GInteractions object
        map <- findOverlaps(GRanges(anchors_to_gene$range), .gi)
        anchors_to_gene <- anchors_to_gene[queryHits(map),]
        anchors_to_gene$which <- subjectHits(map)
        df <- anchors_to_gene[c("which", .cols)]
        df <- chop(df, all_of(.cols))
        mapped_list <- lapply(
            df[.cols], function(x) as(as.list(x), "CompressedList")
        )
        ## Now replicate the data as it would be passed WITH mappings
        .gi <- .gi[df$which]
        mcols(.gi)[.cols] <- mapped_list
    }

    ## If .gi has mappings already (or has them from above)
    mapped_df <- mcols(.gi)[.cols]
    mapped_df[["subjectHits"]] <- seq_along(.gi)
    mapped_df <- as.data.frame(mapped_df)
    mapped_df <- mutate_all(mapped_df, vec_proxy) # Remove AsIs attributes!!!
    mapped_df <- unnest(mapped_df, everything())
    ## Now map back to the original ranges
    gr_to_gi <- as.data.frame(findOverlaps(.gr, .gi, maxgap = .gr2gi, ...))
    range_df <- data.frame(
        range = as.character(.gr), queryHits = seq_along(.gr)
    )
    ## This will drop any additional columns. They can be replaced back in the
    ## parent function
    range_df <- inner_join(
        range_df, gr_to_gi, by = "queryHits", multiple = "all"
    )
    range_df <- inner_join(
        range_df, mapped_df, by = "subjectHits", multiple = "all"
    )
    range_df <- range_df[c("range", .cols)]
    return(range_df)

}

#' @title Map ranges to all genes within a set distance
#' @param .gr The ranges to map onto
#' @param .genes GRanges object containing gene-level information
#' @param .cols The columns from `.genes` to map onto `.gr`
#' @param .within The maximum distance between ranges & genes
#' @param ... Passed to findOverlaps
#' @return A data.frame
#' @importFrom S4Vectors mcols
#' @importFrom vctrs vec_proxy
#' @importFrom dplyr inner_join mutate_all across distinct bind_rows
#' @importFrom tidyr unnest
#' @importFrom tidyselect everything ends_with
#' @import GenomicRanges
#' @keywords internal
.mapWithin <- function(.gr, .genes, .cols, .within, ...) {
    distance <- c() ## R CMD error avoidance
    if (missing(.genes) | length(.gr) == 0) return(NULL)
    .cols <- intersect(.cols, names(mcols(.genes)))
    stopifnot(length(.cols) > 0)
    mapped_df <- mcols(.genes)[.cols]
    mapped_df[["subjectHits"]] <- seq_along(.genes)
    mapped_df <- as.data.frame(mapped_df)
    mapped_df <- mutate_all(mapped_df, vec_proxy) # Remove AsIs attributes!!!
    mapped_df <- unnest(mapped_df, everything())

    ## First find all direct overlaps
    range_to_genes <- as.data.frame(findOverlaps(.gr, .genes, ...))
    ## Now find the nearest within '.within' and collect all unique mappings
    range_to_nearest <- as.data.frame(distanceToNearest(.gr, .genes, ...))
    range_to_nearest <- subset(range_to_nearest, distance <= .within)
    range_to_genes <- bind_rows(range_to_genes, range_to_nearest)
    range_to_genes <- distinct(range_to_genes, across(ends_with("Hits")))

    ## Now define the mappings back to the original set
    range_df <- data.frame(range = as.character(.gr), queryHits = seq_along(.gr))
    range_df <- inner_join(
        range_df, range_to_genes, by = "queryHits", multiple = "all"
    )
    range_df <- inner_join(
        range_df, mapped_df, by = "subjectHits", multiple = "all"
    )
    range_df[c("range", .cols)]
}

.checkMappingArgs <- function(.gr, .genes, .prom, .enh, .gi, .cols, .distArgs) {

    if (missing(.gr)) stop("query ranges must be provided")
    if (!is(.gr, "GRanges"))
        stop("query object must be a GenomicRanges object")
    if (length(.cols) == 0)
        stop("At least one gene/ID column must be specified")
    msg <- NULL
    colInProm <- colInEnh <- colInGI <- colInGenes <- FALSE

    if (!missing(.prom)) {
        if (!is(.prom, "GRanges")) {
            msg <- c(msg, "Promoters must be GRanges\n")
        } else {
            colInProm <- any(.cols %in% .mcolnames(.prom))
        }
    }

    if (!missing(.enh)) {
        if (!is(.enh, "GRanges")) {
            msg <- c(msg, "Enhancers must be GRanges\n")
        } else {
            colInEnh <- any(.cols %in% .mcolnames(.enh))
        }
    }

    if (!missing(.gi)) {
        if (!is(.gi, "GInteractions")) {
            msg <- c(msg, "Interactions must be a GInteractions object\n")
        } else {
            colInGI <- any(.cols %in% .mcolnames(.gi))
        }
    }

    if (!missing(.genes)) {
        ## Genes can be absent if the columns are provided in the other object
        if (!is(.genes, "GRanges")) {
            msg <- c(msg, "Genes must be provided as GRanges\n")
        } else {
            colInGenes <- any(.cols %in% .mcolnames(.genes))
        }
    }

    if (!any(colInProm, colInEnh, colInGI, colInGenes))
        msg <- c(msg, "No requested columns were able to be found\n")

    nonNumeric <- !is.numeric(unlist(.distArgs))
    if (any(nonNumeric))
        msg <- c(msg, "All distance arguments must be numeric\n")

    if (is.null(msg)) return(TRUE)
    message(msg)
    FALSE

}
