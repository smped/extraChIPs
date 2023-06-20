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
#' @usage data(grch37.cytobands)
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

#' @usage data(grch38.cytobands)
#' @rdname cytobands
"grch38.cytobands"


#' @title Datasets for an example region
#'
#' @description Various example datasets for demonstrating analysis and
#' visualisation strategies.
#' Generation of all datasets is documented in
#' `system.file("script/ex_datasets.md", package = "extraChIPs")`
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
#' @usage data(ex_trans)
"ex_trans"

#' @rdname ex_datasets
#' @usage data(ex_genes)
"ex_genes"

#' @rdname ex_datasets
#' @usage data(ex_prom)
"ex_prom"

#' @rdname ex_datasets
#' @usage data(ex_hic)
"ex_hic"

#' @title Datasets for the Fixed-Width Vignette
#'
#' @description GRangesList of peaks and SummarizedExperiment of counts
#' All were saved during initial vignette preparation at
#' https://github.com/smped/extraChIPs_vignette/blob/main/differential_signal_fixed.Rmd
#'
#' @examples
#' data(se)
#' se
#' data(peaks)
#' peaks
#' @name fixed_width_datasets
#' @rdname fixed_width_datasets
#' @usage data(se)
"se"

#' @rdname fixed_width_datasets
#' @usage data(peaks)
"peaks"
