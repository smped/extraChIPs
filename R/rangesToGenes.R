#' @title Map Genomic Ranges to Genes
#'
#' @description Map Genomic Ranges to Genes using Multiple Data Sources
#'
#' @details
#' This function is able to incorporate feature-level information and HiC
#' interactions to enable better mapping of regions to genes. For GRanges, the
#' following sequential strategy is used.
#'
#' \enumerate{
#'   \item Ranges overlapping a promoter are assigned to that gene
#'   \item Ranges overlapping an enhancer are assigned to all genes within
#'   a specified distance
#'   \item Ranges overlapping a HiC interaction are assigned to all genes
#'   connected by the interaction
#'   \item Ranges with no gene assignment from the previous steps are assigned
#'   to the nearest gene within a specified distance
#' }
#'
#' For HiC Interactions, the above strategy is again used, excluding the HiC
#' step and taking each anchor as a distinct range. Mappings for both anchors
#' are then combined for interaction-specific mappings
#'
#' @return
#' A GRanges object with added mcols as specified
#'
#' @param gr GRanges object with ranges to be mapped to genes
#' @param prom GRanges object defining promoters
#' @param enh GRanges object defining Enhancers
#' @param hic GInteractions object defining interactions. Mappings from
#' interactions to genes should be performed as a separate prior step.
#' @param genes GRanges object containing genes (or any other nominal feature)
#' to be assigned
#' @param cols Column names to be assigned as mcols in the output. Columns
#' must be present in both `hic` and `genes` if using both objects as
#' references. Any columns not found in reference objects will be ignored
#' @param maxgap The maximum distance between range and features to be
#' considered as a direct overlap
#' @param maxdist The maximum distance between ranges/enhancers and genes for
#' mapping any otherwise unmapped ranges
#' @param ... Not used
#'
#' @importFrom S4Vectors mcols
rangesToGenes <- function(
  gr, prom, enh, hic, genes, cols = c("gene_id", "gene_name", "symbol"),
  maxgap = 0, maxdist = 1e5, ...
) {
  checkArgs <- .checkMappingArgs(gr, prom, enh, hic, genes, cols)
  stopifnot(checkArgs)

  ## Remove any columns that clash in the original object
  mcols(gr) <- mcols(gr)[setdiff(colnames(mcols(gr)), cols)]

  ## First perform the direct mapping to genes
  gr <- .mapGenesByFeature(gr, prom, enh, genes, cols, maxgap, maxdist, ...)

  ## Now map using HiC data

  ## Map any unmapped genes


}

#' @importFrom IRanges overlapsAny
#' @importFrom dplyr case_when
#' @importFrom GenomicRanges GRanges
.mapGenesByFeature <- function(
  .gr, .prom, .enh, .genes, .cols, .maxgap, .maxdist, ...
){

  ## If no mappings are possible, just return the original object
  if (missing(.genes)) return(.gr)
  if (missing(.prom) & missing(.enh)) return(.gr)
  ## Set up empty ranges if any features are missing
  if (missing(.prom)) .prom <- GRanges()
  if (missing(.enh)) .enh <- GRanges()
  ol <- case_when(
    overlapsAny(.gr, .prom, maxgap = .maxgap) ~ "Promoter",
    overlapsAny(.gr, .enh, maxgap = .maxgap) ~ "Enhancer",
    TRUE ~ "None"
  )
  grl <- as.list(split(.gr, f = ol))
  ## Now map promoters by direct overlap

  ## Map Enhancers to genes with .maxdist

  ## return the list in identical form, but with added columns as requested

}


#' @importFrom S4Vectors mcols
.checkMappingArgs <- function(.gr, .prom, .enh, .hic, .genes, .cols){

  if (missing(.hic) & missing(.genes))
    stop("Either 'hic' or 'genes' must be supplied\n")

  msg <- hiccols <- genecols <- c()
  if (!missing(.hic)) {

    if (!is(.hic, "GInteractions")) {
      msg <- c(msg, "'hic' must be provided as GInteractions\n")
    } else {
      hiccols <- intersect(.cols, colnames(mcols(.hic)))
    }

  }

  if (!missing(.genes)) {

    if (!is(.genes, "GRanges")) {
      msg <- c(msg, "'genes' must be a GRanges object\n")
    } else {
      genecols <- intersect(.cols, colnames(mcols(.genes)))
      if (!is.null(hiccols))
        if (length(genecols) != length(hiccols))
          msg <- c(
            msg,
            "All valid columns in 'cols' must be in both 'hic' and 'genes'\n"
          )
    }

  }

  if (length(c(genecols, hiccols)) == 0)
    msg <- c(msg, "No columns found in 'cols' for obtaining gene information\n")

  if (!is(.gr, "GRanges")) msg <- c(msg, "'gr' must be a GRanges object\n")

  if (!missing(.prom))
    if (!is(.prom, "GRanges"))
      msg <- c(msg, "'prom' must be a GRanges object\n")

  if (!missing(.enh))
    if (!is(.enh, "GRanges")) msg <- c(msg, "'enh' must be a GRanges object\n")

  if (is.null(msg)) return(TRUE)
  message(msg)
  FALSE

}
