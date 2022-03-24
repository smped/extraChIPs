#' @title Cytogenetic bands
#'
#' @description Cytogenetic bands for GRCh37/hg19 and GRCh38/hg38
#'
#' @format Cytogenetic bands for standard chromosomes from GRCh37,in the format
#' required by \link[Gviz]{IdeogramTrack}. A data.frame with 5 columns:
#' \describe{
#'  \item{chrom}{Chromosome}
#'  \item{chromStart}{Starting position for each cytogenetic band}
#'  \item{chromEnd}{End position for each cytogenetic band}
#'  \item{name}{Name for each band, e.g. p.36.33}
#'  \item{gieStain}{Staining pattern}
#' }
#'
#' @examples
#' data(grch37.cytobands)
#' head(grch37.cytobands)
#'
#' data(grch38.cytobands)
#' head(grch38.cytobands)
#'
#' @source \url{https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz}
#' @source \url{https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz}
#' @name cytobands
#' @rdname cytobands
"grch37.cytobands"

#' @rdname cytobands
"grch38.cytobands"


#' @title Datasets for an example region
#'
#' @description Various example datasets for demonstrating visualisation
#' strategies.
#'
#' \describe{
#'   \item{ex_genes}{Simple GRanges object with complete ranges for each gene}
#'   \item{ex_trans}{Exon & transcript level information prepared for plotting
#'   with `Gviz` or `plotHFGC()`}
#'   \item{ex_prom}{Regions defined as promoters}
#'   \item{ex_hic}{Example HiC interactions}
#' }
#'
#' @format GRanges and GInteractions objects
#'
#' All annotations are from GRCh37
#' @examples
#' data(ex_trans)
#' ex_trans
#' @name ex_datasets
#' @rdname ex_datasets
"ex_trans"

#' @rdname ex_datasets
"ex_genes"

#' @rdname ex_datasets
"ex_prom"

#' @rdname ex_datasets
"ex_hic"
