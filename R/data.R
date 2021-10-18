#' @title Cytobands for GRCh37/hg19
#'
#' @description Cytogenetic bands from GRCh37/hg19
#'
#' @format Cytogenetic bands for standard chromosomes from GRCh37,in the format
#' required by \link[Gviz]{IdeogramTrack}. A data.frame with 862 rows and
#' 5 columns:
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
#' @source \url{https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz}
"grch37.cytobands"

#' @title Cytobands for GRCh38/hg38
#'
#' @description Cytogenetic bands from GRCh38/hg38
#'
#' @format Cytogenetic bands for standard chromosomes from GRCh38,in the format
#' required by \link[Gviz]{IdeogramTrack}. A data.frame with 862 rows and
#' 5 columns:
#' \describe{
#'  \item{chrom}{Chromosome}
#'  \item{chromStart}{Starting position for each cytogenetic band}
#'  \item{chromEnd}{End position for each cytogenetic band}
#'  \item{name}{Name for each band, e.g. p.36.33}
#'  \item{gieStain}{Staining pattern}
#' }
#'
#' #' @examples
#' data(grch38.cytobands)
#' head(grch38.cytobands)
#'
#' @source \url{https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz}
"grch38.cytobands"
