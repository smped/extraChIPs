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
#' @examples
#' data(grch38.cytobands)
#' head(grch38.cytobands)
#'
#' @source \url{https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz}
"grch38.cytobands"

#' @title Genes in an example region
#'
#' @description Genes in an example region formatted for Gviz
#'
#' @format Gene structure as a GRanges object, in the correct structure for
#' plotting with plotHFGC.
#' This is a GRangesList with two elements: 'Unchanged' and 'Up' representing
#' genes that were considered as unchanged or upregulated in a separate RNA-seq
#' analysis.
#' Each element is pre-formatted for compatability with the package `Gviz` and
#' has the columns
#' \describe{
#'   \item{type}{The type of feature. All here are exons}
#'   \item{gene}{The Ensembl gene ID}
#'   \item{exon}{An arbitrary exon ID}
#'   \item{transcript}{The Ensembl transcript ID}
#'   \item{symbol}{The gene nme (or symbol)}
#'   }
#'
#' All annotations are from GRCh37
#' @examples
#' data(ex_trans)
#' ex_trans
"ex_trans"

#' @title Basic GRanges containing four genes
#'
#' @description Genes in an example region
#'
#' @format Simple GRanges object to demonstrate mapping from ranges to genes
#'
#' All annotations are from GRCh37
#' @examples
#' data(ex_genes)
#' ex_genes
"ex_genes"

#' @title Basic GRanges containing promoters
#'
#' @description Promoters from an example region
#'
#' @format Simple GRanges object to demonstrate using features in various steps
#' implemented in this package
#'
#' All annotations are from GRCh37
#' @examples
#' data(ex_prom)
#' ex_prom
"ex_prom"

#' @title Basic GInteractions object
#'
#' @description Contains a single interaction fo demonstration purposes
#'
#' @format Simple GInteractions object to demonstrate using features in various
#' steps implemented in this package
#'
#' All annotations are from GRCh37
#' @examples
#' data(ex_hic)
#' ex_hic
"ex_hic"
