#' extraChIPs: A package for enabling and extending ChIP-Seq analysis
#'
#' The package provides three categories of important functions:
#' Range-based, Visualisation and Convenience functions, with most centred
#' around GenomicRanges objects
#'
#' @section Range-based Functions:
#'
#' Many of the range-based functions included in this package have a focus on
#' retaining the `mcols` information whilst manipulating the ranges, such as
#' [reduceMC()] which not only reduces the Ranges, but collapses the `mcols`
#' into vectors or [IRanges::CompressedList] objects.
#' Key function from this group are:
#'
#'  * [reduceMC()], [setdiffMC()], [intersectMC()], [unionMC()], [distinctMC()]
#'  and [chopMC()]
#'  * [bestOverlap()] and [propOverlap()] provide simple output easily able to
#'  be added as a column within the `mcols` element
#'  * [as_tibble()] coerces a GRanges object to a [tibble::tibble].
#'  * [colToRanges()] enables parsing of a single column to a GRanges object,
#'  setting all other columns as the `mcols` element.
#'  * [stitchRanges()] merges nearby ranges setting barrier ranges which cannot
#'  be crossed when merging
#'  * [partitionRanges()] break apart one set of ranges by another
#'  * [dualFilter()] filters ranges from sliding windows using a guide set of
#'  reference ranges where signal is confidently expected
#'  * [mergeByCol()] merges overlapping ranges, as produced by sliding windows
#'  * [mapByFeature()] is able to map a set of GRanges to the most appropriate
#'  gene, using any optional combination of promoters, enhancers and HiC
#'  interactions
#'  * [grlToSE()] takes selected columns from a GRangesList and sets them as
#'  assays within a [SummarizedExperiment::RangedSummarizedExperiment] object.
#'  Used for combining peak intensities or results across multiple ChIP targets.
#'
#' @section Visualisation Functions:
#'
#'  * [plotHFGC()] is a wrapper to Gviz plotting functions, able to take any
#'  combination of HiC, Features, Genes and Coverage (i.e. BigWig) and plot a
#'  specified range.
#'  * [plotOverlaps()] visualises overlapping ranges as an UpSet plot or Venn
#'  Diagram
#'  * [plotProfileHeatmap()] plots the average signal around a set of ranges,
#'   as prepared by [getProfileData()]
#'  * [plotPie()] and [plotSplitDonut()] enable simple comparison across
#'  multiple annotation columns within a GRanges object.
#'  * [plotAssayDensities()], [plotAssayPCA()] and [plotAssayRle()] provide
#'  simple interfaces to plotting key values from a
#'  [SummarizedExperiment::RangedSummarizedExperiment].
#'
#' @section Convenience Functions:
#'   * [fitAssayDiff()] enables differential signal analysis on a
#'   SummarizedExperiment object
#'   * [collapseGenes()] prints a vector of genes for an rmarkdown document,
#'   using italics.
#'   * [importPeaks()] imports large numbers of broadPeak or narrowPeak files
#'   * [makeConsensus()] forms consensus peaks from overlapping ranges within a
#'   GRangesList()
#'   * [voomWeightsFromCPM()] allows creation of an [limma::EList-class] object
#'   as would be created from counts by [limma::voom()], but using
#'   [edgeR::cpm()] values as input.
#'
#' @author
#' Stephen Pederson
#'
#' @docType package
#' @name extraChIPs-package
NULL
