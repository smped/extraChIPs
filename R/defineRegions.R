#' @title Define Genomic Regions Based on Gene Annotations
#'
#' @description
#' Use gene, transcript and exon annotations to define genomic regions
#'
#' @param genes,transcripts,exons GRanges objects defining each level of
#' annotation
#' @param promoter Numeric vector defining upstream and/or downstream distances
#' for a promoter. Passing a single value will define a symmetrical promoter
#' The first value represents the upstream range
#' @param upstream The distance from a TSS defining an upstream promoter
#' @param intron logical(1) Separate gene bodies into introns and exons. If
#' `intron = FALSE` gene bodies will simply be defined as gene bodies
#' @param proximal Distance from a gene to be considered a proximal intergenic
#' region. If set to 0, intergenic regions will simply be considered as
#' uniformly intergenic
#' @param cols Column names to be retained from the supplied annotations
#' @param ... Not used
#'
#' @details
#' Using GRanges annotated as genes, transcripts and exons this function will
#' define ranges uniquely assigned to a region type using a hierarchical
#' process. By default, these region types will be (in order) 1) Promoters, 2)
#' Upstream Promoters, 3) Exons, 4) Introns, 5) Proximal Intergenic and 6)
#' Distal Intergenic.
#'
#' Setting `intron = FALSE` will replace introns and exons with a generic "Gene
#' Body" annotation.
#' Setting `proximal = 0` will return all intergenic regions (not previously
#' annotated as promoters or upstream promoters) to an "Intergenic" category
#'
#' Notably, once a region has been defined, it is excluded from all subsquent
#' candidate regions.
#'
#' Any columns matching the names provided in cols will be returned, and it is
#' assumed that the gene/transcript/exon ranges will contain informative columns
#' in the `mcols()` element.
#'
#' @return A GRangesList
#'
#' @examples
#'
#' ## Define two exons for two transcripts
#' sq <- Seqinfo(seqnames = "chr1", seqlengths = 50000)
#' e <- c("chr1:20001-21000", "chr1:29001-29950", "chr1:22001-23000", "chr1:29001-30000")
#' e <- GRanges(e, seqinfo = sq)
#' mcols(e) <- DataFrame(gene_id = "Gene1", transcript_id = c("Trans1", "Trans1", "Trans2", "Trans2"))
#' mcols(e) <- DataFrame(
#'   gene_id = "Gene1",
#'   transcript_id = c("Trans1", "Trans1", "Trans2", "Trans2")
#' )
#'
#' ## Define the transcript ranges
#' t <- unlist(endoapply(split(e, e$transcript_id), range))
#' t$gene_id <- "Gene1"
#' t$transcript_id <- names(t)
#' names(t) <- NULL
#'
#' ## Summarise to gene level
#' g <- reduceMC(t)
#' g$transcript_id <- NA_character_
#'
#' ## Now annotate the regions
#' regions <- defineRegions(genes = g, transcripts = t, exons = e)
#' sort(unlist(regions))
#'
#' ## Alternatively, collpse gene body and intergenic ranges
#' regions <- defineRegions(
#'   genes = g, transcripts = t, exons = e, intron = FALSE, proximal = 0
#' )
#' sort(unlist(regions))
#'
#'
#' @import GenomicRanges
#' @importFrom S4Vectors mcols 'mcols<-'
#' @export
defineRegions <- function(
    genes, transcripts, exons, promoter = c(2500, 500), upstream = 5000,
    intron = TRUE, proximal = 10000,
    cols = c("gene_id", "gene_name", "transcript_id", "transcript_name"), ...
) {

  ## Sanity checks
  stopifnot(
    all(is(genes, "GRanges"), is(transcripts, "GRanges"), is(exons, "GRanges"))
  )
  sq <- seqinfo(genes)
  stopifnot(is.numeric(c(promoter, upstream, proximal)))
  all_cols <- c(
    colnames(mcols(genes)), colnames(mcols(transcripts)), colnames(mcols(exons))
  )
  stopifnot(all(table(all_cols) == 3))
  cols <- c("region", intersect(cols, all_cols))

  ## Promoters
  promoter <- abs(rep_len(promoter, 2))
  prom <- promoters(transcripts, promoter[1], promoter[2])
  prom$region <- paste0("Promoter (-", promoter[1], "/+", promoter[2], ")")
  mcols(prom) <- mcols(prom)[cols]
  prom <- reduceMC(prom, ignore.strand = TRUE)
  out <- GRangesList(promoter = prom)

  ## Upstream Promoters
  up_prom <- promoters(transcripts, upstream[1], downstream = 0)
  up_prom$region <- paste0("Upstream Promoter (-", upstream / 1e3, "kb)")
  mcols(up_prom) <- mcols(up_prom)[cols]
  up_prom <- setdiffMC(up_prom, prom, ignore.strand = TRUE)
  out$upstream_promoter <- reduceMC(up_prom)

  ## Gene Bodies
  if (intron) {
    ex <- exons
    ex$region <- "Exon"
    mcols(ex) <- mcols(ex)[cols]
    ex <- setdiffMC(ex, unlist(out), ignore.strand = TRUE)
    out$exon <- reduceMC(ex)
    introns <- setdiffMC(genes, unlist(out), ignore.strand = TRUE)
    introns$region <- "Intron"
    mcols(introns) <- mcols(introns)[cols]
    out$intron <- reduceMC(introns)
  } else {
    gn <- setdiffMC(genes, unlist(out), ignore.strand = TRUE)
    gn$region <- "Gene Body"
    mcols(gn) <- mcols(gn)[cols]
    out$gene_body <- reduceMC(gn)
  }

  if (proximal > 0) {
    ## Only define proximal/distal if proximal > 0
    prox <- resize(genes, width = width(genes) + 2 * proximal, fix = 'center')
    prox <- setdiffMC(prox, unlist(out), ignore.strand = TRUE)
    prox$region <- paste0("Intergenic (<", proximal / 1e3, "kb)")
    mcols(prox) <- mcols(prox)[cols]
    out$proximal_intergenic <- reduceMC(prox)
    distal <- setdiffMC(GRanges(sq), unlist(out))
    distal$region <- paste0("Intergenic (>", proximal / 1e3, "kb)")
    out$distal_intergenic <- distal
  } else {
    inter <- setdiffMC(GRanges(sq), unlist(out))
    inter$region <- "Intergenic"
    out$intergenic <- inter
  }

  ## Remove any empty elements
  keep <- vapply(out, length, integer(1)) > 0
  out[keep]

}

