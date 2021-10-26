#' @title Map Genomic Ranges to Genes
#'
#' @description Map Genomic Ranges to Genes using Key Features
#'
#' @details
#' This function is able to incorporate feature-level information and long-range
#' interactions to enable better mapping of regions to genes. For GRanges, the
#' following sequential strategy is used.
#'
#' \enumerate{
#'   \item Ranges overlapping a promoter are assigned to that gene
#'   \item Ranges overlapping an enhancer are assigned to all genes within
#'   a specified distance
#'   \item Ranges overlapping a long-range interaction are assigned to all genes
#'   connected by the interaction
#'   \item Ranges with no gene assignment from the previous steps are assigned
#'   to the nearest gene within a specified distance
#' }
#'
#' For long-range interactions, the above strategy is again used, taking each
#' anchor as a distinct range. Mappings for both anchors are then combined for
#' interaction-specific mappings
#'
#' @return
#' A GRanges object with added mcols as specified
#'
#' @param gr GRanges object with ranges to be mapped to genes
#' @param prom GRanges object defining promoters
#' @param enh GRanges object defining Enhancers
#' @param gi GInteractions object defining interactions. Mappings from
#' interactions to genes should be performed as a separate prior step.
#' @param genes GRanges object containing genes (or any other nominal feature)
#' to be assigned
#' @param cols Column names to be assigned as mcols in the output. Columns
#' must be present in both `gi` and `genes` if using both objects as
#' references. Any columns not found in reference objects will be ignored
#' @param range2prom The maximum permissible distance between a query range and
#' any ranges defined as promoters
#' @param range2enh The maximum permissible distance between a query range and
#' any ranges defined as enhancers
#' @param range2gi The maximum permissible distance between a query range and
#' any ranges defined as GInteraction anchors
#' @param prom2gene The maximum permissible distance between a range provided
#' in `prom` and a gene
#' @param enh2gene The maximum permissible distance between a range provided
#' in `enh` and a gene
#' @param gi2gene The maximum permissible distance between a GInteractions
#' anchor (provided in `gi`) and a gene
#' @param range2gene The maximum permissible distance between a query range and
#' genes (for ranges not otherwise mapped)
#' @param ... Passed to findOverlaps nad overlapsAny internally
#'
#' @importFrom S4Vectors mcols subjectHits
#' @importFrom dplyr bind_rows distinct across
#' @importFrom tidyr chop
#' @importFrom tidyselect all_of
#' @importFrom GenomicRanges findOverlaps GRanges
#' @importFrom methods as
#' @importFrom stats setNames
mapByFeature <- function(
  gr, genes, prom, enh, gi, cols = c("gene_id", "gene_name", "symbol"),
  range2prom = 0, range2enh = 0, range2gi = 0, range2gene = 1e5,
  prom2gene = 0, enh2gene = 1e5, gi2gene = 0, ...
) {

  ## Perform mapping in a sequential manner as described in the algorithm
  ## First identify the mappings from gr to features
  ol_prom <- ol_enh <- ol_gi <- rep(FALSE, length(gr))
  if (!missing(prom)) ol_prom <- overlapsAny(gr, prom, maxgap = range2prom, ...)
  if (!missing(enh)) ol_enh <- overlapsAny(gr, enh, maxgap = range2enh, ...)
  if (!missing(gi)) ol_gi <- overlapsAny(gr, gi, maxgap = range2gi, ...)
  ol_none <- !ol_prom & !ol_enh & !ol_gi

  ## 1. Map promoters
  prom_tbl <- .mapFeatures(
    gr[ol_prom,], prom, genes, cols, range2prom, prom2gene, ...
  )
  ## 2. Map enhancers
  enh_tbl <- .mapFeatures(
    gr[ol_enh], enh, genes, cols, range2enh, enh2gene, ...
  )
  ## 3. Add gi mappings
  gi_tbl <- .mapGi(gr[ol_gi], gi, genes, cols,  range2gi, gi2gene, ...)
  ## 4. Map any unmapped regions
  misc_tbl <- .mapWithin(gr[ol_none], genes, cols, range2gene, ...)

  ## Return the complete set of mappings & keep them as unique mappings
  mapped_tbl <- bind_rows(list(prom_tbl, enh_tbl, bind_rows(gi_tbl), misc_tbl))
  cols <- intersect(cols, names(mapped_tbl))
  mapped_tbl <- distinct(mapped_tbl, across(all_of(cols)), .keep_all = TRUE)
  mapped_tbl <- chop(mapped_tbl, all_of(cols))
  ## Prepare for returning as an addition to the original GRanges
  index <- findOverlaps(gr, GRanges(mapped_tbl[["range"]]), ...)
  mapped_tbl <- mapped_tbl[subjectHits(index),]
  mapped_list <- as.list(mapped_tbl[cols])
  list_cols <- vapply(mapped_list, is.list, logical(1))
  mapped_list[list_cols] <- lapply(
    ## The new vctrs-based types are problematic still. coerce via as.list()
    mapped_list[list_cols], function(x) as(as.list(x), "List")
  )
  grl <- as.list(split(gr, f = seq_along(gr) %in% queryHits(index)))
  mcols(grl[["TRUE"]])[cols] <- mapped_list
  empty_list <- vector("list", length(grl[["FALSE"]]))
  unmapped_list <- lapply(mapped_list, function(x) as(empty_list, is(x)[[1]]))
  mcols(grl[["FALSE"]])[cols] <- unmapped_list
  unsorted <- setNames(unlist(GRangesList(grl)), NULL)
  o <- subjectHits(findOverlaps(gr, unsorted))
  stopifnot(length(o) == length(gr))
  unsorted[o,]

}

#' @title Map ranges to genes using features as an anchor
#' @param .gr The ranges to map onto
#' @param .feat Features to use for mapping
#' @param .genes GRanges object containing gene-level information
#' @param .cols The columns from `.genes` to map onto `.gr`
#' @param .gr2feat The maximum distance between ranges and features
#' @param .feat2gene The maximum distance between features & genes
#' @param ... Passed to findOverlaps
#' @return A data.frame
#' @importFrom S4Vectors mcols
#' @importFrom GenomicRanges findOverlaps
#' @importFrom dplyr inner_join left_join mutate_all
#' @importFrom vctrs vec_proxy
#' @importFrom tidyr unnest
#' @importFrom tidyselect everything
.mapFeatures <- function(.gr, .feat, .genes, .cols, .gr2feat, .feat2gene, ...) {

  if (missing(.feat)) return(NULL)
  if (length(.feat) == 0 | length(.gr) == 0) return(NULL)
  if (!missing(.genes)) {
    .cols <- intersect(.cols, names(mcols(.genes)))
  } else {
    .cols <- intersect(.cols, names(mcols(.feat)))
  }
  stopifnot(length(.cols) > 0)
  if (!all(.cols %in% names(mcols(.feat)))) {
    ## Map to genes if any requested column is absent from any existing mapping
    feat_to_gene <- findOverlaps(.feat, .genes, maxgap = .feat2gene, ...)
    feat_to_gene <- as.data.frame(feat_to_gene)
    feat_df <- data.frame(queryHits = seq_along(.feat))
    feat_df <- left_join(feat_df, feat_to_gene, by = "queryHits")
    .genes$subjectHits <- seq_along(.genes)
    gene_df <- as.data.frame(mcols(.genes)[c("subjectHits", .cols)])
    feat_df <- left_join(feat_df, gene_df, by = "subjectHits")
    feat_df <- feat_df[c("queryHits", .cols)]
    names(feat_df) <- gsub("queryHits", "subjectHits", names(feat_df))
  } else {
    ## Setup as the same format otherwise
    feat_df <- mcols(.feat)[.cols]
    feat_df[["subjectHits"]] <- seq_along(.feat)
    feat_df <- as.data.frame(feat_df)
    feat_df <- mutate_all(feat_df, vec_proxy) # Remove an AsIs attributes!!!
    feat_df <- unnest(feat_df, everything())
  }

  ## Map ranges to featoters
  range_to_feat <- findOverlaps(.gr, .feat, maxgap = .gr2feat, ...)
  range_to_feat <- as.data.frame(range_to_feat)
  range_df <- as_tibble(.gr)
  range_df$queryHits <- seq_along(.gr)
  ## This will drop any additional columns. They can be replaced back in the
  ## parent function
  range_df <- inner_join(range_df, range_to_feat, by = "queryHits")
  range_df <- left_join(range_df, feat_df, by = "subjectHits")
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
#' @return List (length 2) of tibbles
#' @importFrom S4Vectors mcols
#' @importFrom vctrs vec_proxy
#' @importFrom GenomicRanges findOverlaps
#' @importFrom InteractionSet anchors
#' @importFrom dplyr inner_join left_join mutate_all
#' @importFrom tidyr unnest
#' @importFrom tidyselect everything
.mapGi <- function(.gr, .gi, .genes, .cols, .gr2gi, .gi2gene , ...) {

  if (missing(.gi) | length(.gr) == 0) return(NULL)
  if (!missing(.genes)) {
    .cols <- intersect(.cols, names(mcols(.genes)))
  } else {
    .cols <- intersect(.cols, names(mcols(.gi)))
  }
  stopifnot(length(.cols) > 0)

  if (!all(.cols %in% names(mcols(.gi)))) {
    .gi <- subsetByOverlaps(.gi, .gr)
    anchors_to_gene <- lapply(
      anchors(.gi), .mapWithin, .genes = .genes, .cols = .cols, .within = 0
    ) # This returns the anchor range, not .gr
    anchors_to_gene <- bind_rows(anchors_to_gene)
    anchors_to_gene <- mutate_all(anchors_to_gene, vec_proxy)
    anchors_to_gene <- unnest(anchors_to_gene, everything())
    anchors_to_gene$which <- subjectHits(findOverlaps(GRanges(anchors_to_gene$range), .gi))
    tbl <- left_join(tibble(which = seq_along(.gi)), anchors_to_gene)[c("which", .cols)]
    tbl <- chop(tbl, all_of(.cols))
    stopifnot(length(.gi) == nrow(tbl))
    mapped_list <- lapply(tbl[.cols], function(x) as(as.list(x), "List"))
    mcols(.gi)[.cols] <- mapped_list
  }

  ## If we have all of the mappings already, just use them
  mapped_df <- mcols(.gi)[.cols]
  mapped_df[["subjectHits"]] <- seq_along(.gi)
  mapped_df <- as.data.frame(mapped_df)
  mapped_df <- mutate_all(mapped_df, vec_proxy) # Remove AsIs attributes!!!
  mapped_df <- unnest(mapped_df, everything())
  range_to_gi <- as.data.frame(findOverlaps(.gr, .gi, maxgap = .gr2gi, ...))
  range_df <- as_tibble(.gr)
  range_df$queryHits <- seq_along(.gr)
  ## This will drop any additional columns. They can be replaced back in the
  ## parent function
  range_df <- inner_join(range_df, range_to_gi, by = "queryHits")
  range_df <- left_join(range_df, mapped_df, by = "subjectHits")
  range_df <- range_df[c("range", .cols)]
  return(range_df)
  ## Looks right. Needs more testing

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
#' @importFrom GenomicRanges findOverlaps
#' @importFrom dplyr inner_join left_join mutate_all
#' @importFrom tidyr unnest
#' @importFrom tidyselect everything
.mapWithin <- function(.gr, .genes, .cols, .within, ...) {
  if (missing(.genes) | length(.gr) == 0) return(NULL)
  .cols <- intersect(.cols, names(mcols(.genes)))
  stopifnot(length(.cols) > 0)
  mapped_df <- mcols(.genes)[.cols]
  mapped_df[["subjectHits"]] <- seq_along(.genes)
  mapped_df <- as.data.frame(mapped_df)
  mapped_df <- mutate_all(mapped_df, vec_proxy) # Remove AsIs attributes!!!
  mapped_df <- unnest(mapped_df, everything())
  range_to_genes <- findOverlaps(.gr, .genes, maxgap = .within, ...)
  range_to_genes <- as.data.frame(range_to_genes)
  range_df <- as_tibble(.gr)
  range_df$queryHits <- seq_along(.gr)
  ## This will drop any additional columns. They can be replaced in the
  ## parent function
  range_df <- inner_join(range_df, range_to_genes, by = "queryHits")
  range_df <- left_join(range_df, mapped_df, by = "subjectHits")
  range_df[c("range", .cols)]
}
