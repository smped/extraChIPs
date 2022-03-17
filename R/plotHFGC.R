#' @title Plot a Genomic Region showing HiC, Features, Genes and Coverage
#'
#' @description Plot a region with showing HiC, Features, Genes and Coverage
#'
#' @details
#' Convenience function for plotting a common set of tracks. All tracks are
#' optional. For more fine control, users are advised to simply use Gviz
#' directly.
#'
#' The primary tracks defined in this function are H (HiC), F (features), G
#' (genes), and C (coverage). Axis and Ideogram tracks are an additional part of
#' this visualisation
#'
#' Use all tracks specific to this dataset to generate a simple visualisation.
#' In descending order the tracks displayed will be:
#' \enumerate{
#'   \item HiC Interactions (if supplied)
#'   \item Regulatory features
#'   \item Genes/genes
#'   \item Coverage tracks as supplied
#' }
#'
#' All tracks are optional and will simply be omitted if no data is supplied.
#'
#' If wanting a single track of genes, simply pass a GRanges object in the
#' format specified for a \link[Gviz]{GeneRegionTrack}. Passing a GRangesList
#' with the same format will yield an individual track for each list element,
#' with each track shown by default as a separate colour. This can be used for
#' showing Up/Down-regulated genes, or Detected/Undetected genes.
#'
#' If passing a BigWigFileList for the coverage track, each file within the
#' object will be drawn on a separate track. If passing a list of BigWigFileList
#' objects, each list element will be drawn as a single track with the
#' individual files within each BigWigFileList overlaid within each track.
#'
#' Cytogenetic band information must be in the structure required by
#' \link[Gviz]{IdeogramTrack}, with data for both GRCh37 and GRCh38 provided in
#' this package (\link{grch37.cytobands}, \link{grch38.cytobands}).
#'
#' A highlight overlay over the GRanges provided as the `gr` argument will be
#' added if a colour is provided. If set to NULL, no highlight will be added.
#'
#' @examples
#' gr <- GRanges("chr1:11869-12227")
#' feat_gr <- GRangesList(
#'   Promoter = GRanges("chr1:11800-12000"),
#'   Enhancer = GRanges("chr1:13000-13200")
#' )
#' hic <- InteractionSet::GInteractions(feat_gr$Promoter, feat_gr$Enhancer)
#' genes <- c("chr1:11869-12227:+", "chr1:12613-12721:+", "chr1:13221-14409:+")
#' genes <- GRanges(genes)
#' mcols(genes) <- DataFrame(
#'   feature = "exon", gene = "ENSG00000223972", exon = 1:3,
#'   transcript = "ENST00000456328", symbol = "DDX11L1"
#' )
#' data(grch37.cytobands)
#' plotHFGC(
#'   gr, hic = hic, features = feat_gr, genes = genes,
#'   zoom = 2, cytobands = grch37.cytobands, rotation.title = 90,
#'   featurecol = c(Promoter = "red", Enhancer = "yellow")
#' )
#'
#' @return
#' A Gviz object
#'
#' @param gr The range(s) of interest. Must be on a single chromosome
#' @param hic Any HiC interactions to be included as a GenomicInteractions
#' object. If not supplied, no HiC track will be drawn.
#' @param features A named GRangesList object containing regulatory features in
#' each list element. Features will be drawn on a single track with colours
#' matching those provided in `featurecol`. If not included, no feature track
#' will be drawn
#' @param genes A GRanges object with exon structure for each transcript/gene.
#' If not included, no track will be drawn for gene/transcript structure
#' @param coverage A named list of BigWigFileList objects containing the
#' primary tracks to show coverage for. Each list element will be drawn on a
#' separate track, with elements within each BigWigFileList shown on the same
#' track. List names will become track names. Alternatively, a single
#' BigWigFileList will plot all individual files on separate tracks.
#' If not included, no coverage tracks will be drawn.
#' @param annotation Annotations for the coverage track(s).
#' A single GRangesList if coverage is a BigWigListList.
#' If coverage is supplied as a list of BigWigFileLists, a named list of
#' GRangesList objects for each coverage track being annotatated. Names must
#' match those given for coverage.
#' @param zoom Multiplicative factor for zooming in and out
#' @param shift Shift the plot. Applied after zooming
#' @param axistrack logical. Add an AxisTrack()
#' @param cytobands Cytogenetic bands to be displayed on each chromosome
#' @param max The maximum width of the plotting region. Given that the width of
#' the final plotting window will be determined by any HiC interactions, this
#' argument excludes any interactions beyond this distance. Plotting can be
#' somewhat slow if any long range interactions are included. Ignored if no HiC
#' interactions are supplied.
#' @param fontsize Applied across all tracks
#' @param hiccol list with names `"anchors"` and `"interactions"`. Colours
#' are passed to these elements
#' @param featurecol Named vector (or list) of colours for each feature
#' @param genecol Named vector (or list) of colours for each gene category
#' @param annotcol Colours matching the coverage annotations
#' @param coverage_type The plot type for coverage. Currently only lines ("l")
#' and heatmaps ("heatmap") are supported
#' @param linecol If passing a BigWigFileList to coverage, a vector of colours.
#' If passing a list of BigWigFileList objects to coverage, a list of colours
#' with structure that matches the object being passed to coverage, i.e. a
#' named list of the same length, with elements who's length matches each
#' BigWigFileList. Only used if coverage_type = "l".
#' @param gradient Colour gradient for heatmaps
#' @param highlight Outline colour for the highlight track. Setting this to
#' `NULL` will remove the highlight
#' @param hicsize,featsize,genesize,covsize,annotsize
#' Relative sizes for each track (hic, features, genes, coverage & annotation)
#' @param ylim If a numeric vector, this will be passed to all coverage tracks.
#' Alternatively, a named list of y-limits for each coverage track with names
#' that match those in each element of the coverage list.
#' @param ... Passed to \link[Gviz]{DataTrack} for the **coverage tracks** only.
#' Useful arguments may be things like `legend`
#' @param cex.title Passed to all tracks
#' @param rotation.title Passed to all tracks
#' @param collapseTranscripts Passed to \link[Gviz]{GeneRegionTrack} for the
#' genes track
#'
#' @importFrom methods slot
#' @importFrom InteractionSet anchors
#' @importFrom GenomicInteractions calculateDistances
#' @importFrom GenomicRanges GRangesList resize shift
#' @importFrom grDevices hcl.colors
#' @importFrom Gviz HighlightTrack GenomeAxisTrack plotTracks
#' @importFrom IRanges start end width
#'
#' @export
plotHFGC <- function(
  gr, hic, features, genes, coverage, annotation,
  zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, cytobands,
  coverage_type = c("l", "heatmap"),
  linecol = c(), gradient = hcl.colors(101, "viridis"),
  hiccol = list(anchors = "lightblue", interactions = "red"),
  featurecol, genecol, annotcol, highlight = "blue",
  hicsize = 1, featsize = 1, genesize = 1, covsize = 4, annotsize = 0.5,
  ylim = NULL, ...,
  fontsize = 12, cex.title = 0.8, rotation.title = 0,
  collapseTranscripts = "meta"
) {

  ## Argument checks
  coverage_type <- match.arg(coverage_type)
  checkArgs <- .checkHFGCArgs(
    gr = gr, zoom = zoom, shift = shift, hic = hic, features = features,
    genes = genes, coverage = coverage, annotation = annotation,
    axistrack = axistrack, cytobands = cytobands, max = max, hiccol = hiccol,
    linecol = linecol, genecol = genecol, featurecol = featurecol,
    annotcol = annotcol, type = coverage_type, ylim = ylim
  )
  stopifnot(checkArgs)

  ## Form the HiC track, including all interactions beyond the max
  hic_track <- .makeHiCTrack(
    hic, gr, fontsize, hicsize, cex.title, rotation.title, hiccol
  )

  ## If interactions were found, reset the plot range. This should be the
  ## maximum of all interactions < max or the initial range
  plot_range <- range(gr)
  if (!is.null(hic_track)) {
    hic <- slot(hic_track, "giobject")
    anchors <- anchors(hic[calculateDistances(hic) < max])
    plot_range <- range(unlist(GRangesList(anchors)))
  }
  ## Now resize/shift as required
  plot_range <- resize(plot_range, zoom*width(plot_range), fix = "center")
  plot_range <- shift(plot_range, shift)

  ## Form the IdeogramTrack
  ideo_track <- .makeIdeoTrack(plot_range, cytobands, fontsize)

  ## Form the features track
  feature_track <- .makeFeatureTrack(
    features, plot_range, fontsize, featurecol, featsize, cex.title,
    rotation.title
  )

  ## Form the genes tracks. NB: This will be a list of tracks
  gene_tracks <- .makeGeneTracks(
    genes, plot_range, collapseTranscripts, fontsize, genecol, genesize,
    cex.title, rotation.title
  )

  ## The coverage tracks NB: This will be list of tracks
  cov_tracks <- .makeCoverageTracks(
    coverage, plot_range, fontsize, coverage_type, linecol, gradient, covsize,
    cex.title, rotation.title, ylim, ...
  )
  cov_tracks <- .addAnnotations(
    annotation, plot_range, cov_tracks, coverage, annotcol, annotsize
  )

  ## Add the highlight track if wanted. Include features, genes & coverage
  hl_track <- c(feature_track, gene_tracks, cov_tracks)
  hl_track <- hl_track[!vapply(hl_track, is.null, logical(1))]
  if (!is.null(highlight))
    hl_track <- HighlightTrack(
      trackList = hl_track, range = gr, col = highlight, fill = "#FFFFFF00",
      inBackground = FALSE
    )

  plot_list <- list(ideo_track)
  if (axistrack)
    plot_list <- c(plot_list, GenomeAxisTrack(plot_range, fontsize = fontsize))

  plot_list <- c(plot_list, hic_track, hl_track)
  plotTracks(plot_list, from = start(plot_range), end(plot_range))

}


#' @importFrom Gviz IdeogramTrack
#' @importFrom GenomeInfoDb genome seqnames
#' @importFrom rtracklayer ucscGenomes
.makeIdeoTrack <- function(.gr, .bands, .fontsize) {
  ## Checks have been done in .checkHFGCArgs
  ## Any missing cytoband information will have to generated automatically
  ## This can be either downloaded, or generated as a hack from the gr
  ## If a UCSC genome is provided without cytogenetic information: download
  ## Obviously, this is really slow!!!
  ## Otherwise, just use the seqinfo object
  gen <- as.character(unique(genome(.gr)))
  chr <- as.character(unique(seqnames(.gr)))
  if (missing(.bands)) {
    ## First map any GRCh IDs to UCSC
    grc2ucsc <- c(
      GRCh37 = "hg19", GRCh38 = "hg38", GRCm38 = "mm10", GRCm39 = "mm39",
      GRCz10 = "danRer10", GRCz11 = "danRer11", "Rnor_6.0" = "rn6",
      "mRatBN7.2" = "rn7", "Gallus_gallus-5.0" = "galGal5", GRCg6a = "galGal6"
    )
    if (gen %in% names(grc2ucsc)) gen <- grc2ucsc[gen]
    ## Now check for a valid UCSC name & set to NULL for an automatic download
    avail <- ucscGenomes()[,"db"]
    if (gen %in% avail) {
      .bands <- NULL
    } else {
      ## If no valid genome is available return a NULL
      message("Could not find cytogenetic bands for genome ", gen)
      return(NULL)
    }
  }
  IdeogramTrack(
    chromosome = chr, genome = gen, name = chr, bands = .bands,
    fontsize = .fontsize
  )
}

#' @importFrom GenomicInteractions anchorOne anchorTwo
#' @importFrom GenomicInteractions GenomicInteractions InteractionTrack
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges subsetByOverlaps
.makeHiCTrack <- function(.hic, .gr, .fontsize, .tracksize, .cex, .rot, .col) {

  if (missing(.hic)) return(NULL)

  ## Checks have been performed previously on any provided object & params

  ## Just keep cis interactions on the chromosome of interest
  .hic <- subsetByOverlaps(.hic, .gr)
  if (length(.hic) == 0) return(NULL)
  chr <- as.character(seqnames(.gr))[[1]]
  cis <- seqnames(anchorOne(.hic)) == chr & seqnames(anchorTwo(.hic)) == chr
  if (sum(cis) == 0) return(NULL)
  .hic <- .hic[cis]

  track <- InteractionTrack(
    x = GenomicInteractions(.hic),
    chromosome = chr,
    name = "HiC"
  )
  track@dp@pars$size <- .tracksize
  track@dp@pars$fontsize <- .fontsize
  track@dp@pars$cex.title <- .cex
  track@dp@pars$rotation.title <- .rot
  track@dp@pars$col.anchors.line <- .col[["anchors"]]
  track@dp@pars$col.anchors.fill <- .col[["anchors"]]
  track@dp@pars$col.interactions <- .col[["interactions"]]
  track
}

#' @importFrom methods is
.addAnnotations <- function(
  .annotation, .gr, .cov_tracks, .coverage, .fill, .size
  ) {
  if (missing(.annotation) | is.null(.cov_tracks)) return(.cov_tracks)

  if (is(.coverage, "BigWigFileList")) {
    ## Here we just add another feature track. Input will be a single GRList
    fontsize <- .cov_tracks[[1]]@dp@pars$fontsize
    cex <- .cov_tracks[[1]]@dp@pars$cex.title
    annot_track <- .makeFeatureTrack(
      .annotation, .gr, fontsize, .fill, .size, cex, 0
    )
    annot_track@name <- ""
    tracks <- c(list(annot_track), .cov_tracks)
  }

  if (is(.coverage, "list")) {
    ## Here we just add a feature track before any coverage tracks
    ## Find the common names
    # browser()
    cov2add <- intersect(names(.annotation), names(.coverage))
    if (length(cov2add) == 0) return(.cov_tracks)
    tracks <- lapply(
      .cov_tracks,
      function(x) {
        nm <- x@name
        if (!nm %in% cov2add) return(x)
        fontsize <- x@dp@pars$fontsize
        cex <- x@dp@pars$cex.title
        annot_track <- .makeFeatureTrack(
          .annotation[[nm]], .gr, fontsize, .fill, .size, cex, 0
        )
        annot_track@name <- nm
        c(list(annot_track), x)
      }
    )
  }

  unlist(tracks)

}

#' @importFrom GenomicRanges granges
#' @importFrom IRanges subsetByOverlaps
#' @importFrom Gviz AnnotationTrack
.makeFeatureTrack <- function(
  .features, .gr, .fontsize, .fill, .tracksize, .cex, .rot
) {

  if (missing(.features)) return(NULL)

  if (missing(.fill)) {
    pos <- ((seq_along(.features) - 1) %% 8) + 1 # Use modulo for recursion
    .fill <- brewer.pal(8, "Set2")[pos]
    names(.fill) <- names(.features)
  }
  .fill <- unlist(.fill)[names(.features)]

  .features <- granges(unlist(.features))
  .features$feature <- names(.features)
  .features <- subsetByOverlaps(.features, .gr)
  if (length(.features) == 0) return(NULL)
  AnnotationTrack(
    ## Change the name later
    range = .features, name = "Features", stacking = "full",
    col = "transparent", fill = .fill[.features$feature],
    feature = .features$feature,
    ## Tidy setting this up later & also tidy the colour setting
    ## These can likely go in the @dp@pars by just adding them as a list item
    fontsize = .fontsize, cex.title = .cex,
    size = .tracksize,
    rotation.title = .rot
  )

}

#' @importFrom IRanges subsetByOverlaps
#' @importFrom Gviz GeneRegionTrack
#' @importFrom stringr str_to_title
#' @importFrom RColorBrewer brewer.pal
#' @importFrom methods is
.makeGeneTracks <- function(
  .genes, .gr, .collapse, .fontsize, .col, .tracksize, .cex, .rot
) {

  if (missing(.genes)) return(list(NULL))
  gene <- c() # Avoid an R CMD check error

  if (is(.genes, "GRanges")) {
    ## If a single GRanges object, just return a single track
    if (missing(.col)) .col <- "#FFD58A"
    ids <- subsetByOverlaps(.genes, .gr)$gene
    if (length(ids) == 0) return(list(NULL))

    trackList <- GeneRegionTrack(
      subset(.genes, gene %in% ids), name = "Genes",
      transcriptAnnotation = "symbol", collapseTranscripts = .collapse,
      size = .tracksize, fontsize = .fontsize,
      col = "transparent", fill = .col[[1]],
      cex.title = .cex, rotation.title = .rot
    )
    return(trackList)
  }

  ## Below will only execute if .genes is a GRangesList
  ids <- lapply(.genes, subsetByOverlaps, .gr)
  ids <- unlist( lapply(ids, function(x) x$gene) )
  if (length(ids) == 0) return(list(NULL))
  .genes <- lapply(.genes, subset, gene %in% ids)
  .genes <- .genes[vapply(.genes, length, integer(1)) > 0]

  ## Make a default set of colours for the gene tracks
  if (missing(.col)) {
    n <- 9L # Set1 contains 9 colours
    pos <- ((seq_along(.genes) - 1) %% n) + 1 # Use modulo for recursion
    .col <- structure(brewer.pal(n, "Set1")[pos], names = names(.genes))
  }
  .col <- .col[names(.genes)]

  trackList <- lapply(names(.genes),
    function(x) {
      GeneRegionTrack(
        .genes[[x]], name = str_to_title(x),
        transcriptAnnotation = "symbol",
        collapseTranscripts = .collapse,
        size = .tracksize, fontsize = .fontsize, col = .col[[x]],
        fill = .col[[x]], cex.title = .cex, rotation.title = .rot
      )
    })

  return(trackList)

}

#' @importFrom Gviz DataTrack
#' @importFrom rtracklayer import.bw
#' @importFrom S4Vectors mcols
#' @importFrom GenomicRanges GRangesList
.makeCoverageTracks <- function(
  .coverage, .gr, .fontsize, .type = c("l", "heatmap"), .linecol,
  .gradient, .tracksize, .cex, .rot, .ylim,  ...
) {

  ## Always returns a list
  if (missing(.coverage)) return(list(NULL))
  if (is.null(.coverage)) return(list(NULL))
  .type <- match.arg(.type)

  ## Split a single BigWigFileList. This will give a list of BigWigFileLists
  ## (which is a little unexpected). After this step, we will always(!) have a
  ## named S3 list of BigWigFileList objects
  is_single_list <- is(.coverage, "BigWigFileList")
  if (is_single_list) .coverage <- split(.coverage, f = names(.coverage))
  ## Ensure .linecol is a list of the same structure as .coverage
  .linecol <- .assignColours(.coverage, .linecol, .type)
  ## Ensure the .ylim is a list of the same structure as .coverage
  .ylim <- .assignLimits(.coverage, .ylim)

  tracks <- lapply(
    names(.coverage),
    function(x) {
      ## groups need to be unset if drawing a heatmap
      nm <- names(.coverage[[x]])
      grp <- NULL
      if (.type == "l") grp <- factor(nm, levels = nm)
      ## Now import each file within each list element
      cov <- lapply(nm, function(f) {
        data <- import.bw(.coverage[[x]][[f]], which = .gr)
        colnames(mcols(data)) <- f
        data
      })
      cov <- cov[vapply(cov, length, integer(1)) > 0]
      if (length(cov) == 0) return(NULL)
      DataTrack(
        range = unlist(GRangesList(cov)),
        name = x, groups = grp, type = .type,
        fontsize = .fontsize, col = .linecol[[x]], gradient = .gradient,
        showSampleNames = .type == "heatmap", # Always set for heatmaps
        size = .tracksize, cex.title = .cex, rotation.title = .rot,
        ylim = .ylim[[x]], ...
      )
    }
  )
  tracks
}

.assignLimits <- function(.coverage, .ylim) {

  lim <- list() # An empty list
  if (is.list(.ylim)) {
    ## If .ylim is already a list, just make sure the names are present
    lim <- lapply(.ylim, range)
    if (is.null(names(.ylim))) names(lim) <- names(.coverage)
  }

  ## If .ylim is a simple vector, return a list with this as each element
  if (is.numeric(.ylim)) lim <- lapply(.coverage, function(x) range(.ylim))

  ## If no limits are provided, return a list which matches .coverage
  if (is.null(.ylim)) {
    lim <- vector("list", length(.coverage))
    names(lim) <- names(.coverage)
  }
  lim
}

.assignColours <- function(.coverage, .linecol, .type) {

  ## Coverage will always be a list of BigWigFileLists
  ## If the initial object was a single BWFL, then all will internally
  ## have a length of 1. If this is the case, the required structure is a list
  ## of length == length(coverage), with all elements as a single value
  ## Alternatively, if the original object was a list of BWFL objects, we need
  ## a list of colours with the same names as coverage, and with matching
  ## lengths. This gives the structure in essence:
  ## We need a list with the same names as coverage, and elements which have
  ## matching lengths. Even in the case of a heatmap. However, in this case
  ## all elements will be null
  if (.type == "heatmap") {
    cols <- vector("list", length(.coverage))
    names(cols) <- names(.coverage)
    return(cols)
  }
  if (is.null(.linecol)) {
    ## Assign default colours from GViz using a recursion strategy
    defCols <- c(
      "#0080ff", "#ff00ff", "darkgreen", "#ff0000", "orange", "#00ff00", "brown"
    )
    cols <- lapply(
      .coverage, function(x) {
        n <- length(x)
        ind <- ((seq_len(n) - 1) %% 6) + 1 # Use modulo for recursion if n > 6
        defCols[ind]
      })
    return(cols)
  }
  ## If the function has progressed this far, we are using the supplied colours
  ## If the initial input was a BWFL, a vector will have been provided. Split
  ## this into a named list. If a named list was provided, return the same
  all_single <- all(vapply(.coverage, length, integer(1)) == 1)
  if (all_single) {
    cols <- unlist(.linecol) # In case we have a list
    ## Names will have been checked if they are present
    if (is.null(names(cols))) names(cols) <- names(.coverage)
    return(as.list(cols))
  }
  ## If the coverage was a list of BWFL objects, .linecol must have been a
  ## list which matches the structure exactly. This will have been checked in
  ## the first steps of the function and can be returned as is
  .linecol

}



#' @importFrom methods is
#' @importFrom S4Vectors mcols
#' @importFrom GenomeInfoDb seqnames
.checkHFGCArgs <- function(
  gr, zoom, shift, hic, features, genes, coverage, annotation, axistrack,
  cytobands, max, hiccol, linecol, genecol, featurecol, annotcol, type, ylim
) {

  msg <- c()
  if (!is(gr, "GRanges")) {
    msg <- c(msg, "'gr' must be a GRanges object.\n")
  } else {
    if (length(gr) == 0) msg <- c(msg, "Cannot be an empty range\n")

    seq_names <- as.character(unique(seqnames(gr)))
    if (length(seq_names) > 1)
      msg <- c(msg, "All ranges must be on the same chromosome\n")
  }

  nums <- as.numeric(c(zoom, shift, max))
  if (any(is.na(nums)))
    msg <- c( msg, "zoom/shift/max must be numeric or coercible to numeric\n")

  if (!missing(hic)) {
    if (!is(hic, "GInteractions"))
      msg <- c(msg, "'hic' must be provided as a GInteractions object\n")
    colargs <- c("anchors", "interactions")
    if (!is(hiccol, "list") | !all(colargs %in% names(hiccol)))
      msg <- c(
        msg, paste("'hiccols' must be a list with name:", colargs, "\n")
      )
  }

  if (!missing(features)) {
    if (!is(features, "GRangesList"))
      msg <- c(msg, "'features' must be provided as a GRangesList\n")
    if ("" %in% names(features))
      msg <- c(msg, "All elements of 'features' must be explicitly named\n")
    if (!missing(featurecol)) {
      if (!all(names(features) %in% names(featurecol)))
        msg <- c(msg, "All elements of 'features' must be in 'featurecol'\n")
    }

  }

  if (!missing(genes)) {
    if (is(genes, "GRangesList")) {
      if (!missing(genecol))
        if (!all(names(genes) %in% names(genecol)))
          msg <- c(
            msg,
            "All elements of 'genes' must have a colour named in 'genecol'\n"
          )
      genes <- unlist(genes)
    }
    if (!is(genes, "GRanges")) {
      msg <- c(msg, "genes must be a 'GRanges' or GRangesList object\n")
    } else {
      gene_mcols <- c("gene", "exon", "transcript", "symbol")
      if (!all(gene_mcols %in% colnames(mcols(genes))))
        msg <- c(
          msg, paste(
            "'genes' must have an mcols component with the columns:",
            paste(gene_mcols, collapse = ", "), "\n")
        )
    }
  }

  msg <- .checkCoverage(
    msg, coverage, linecol, type, annotation, annotcol, ylim
  )

  bools <- c(axistrack)
  if (!is.logical(bools)) msg <- c(msg, "axistrack must be logical values")

  if (!missing(cytobands)) {
    band_cols <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")
    if (!is(cytobands, "data.frame"))
      msg <- c(msg, "cytobands must be supplied as a data.frame\n")
    if (!all(band_cols %in% colnames(cytobands)))
      msg <- c(
        msg, paste(
          "Columns in cytobands must be exactly:",
          paste(band_cols, collapse = ", "), "\n")
      )
    if (!seq_names %in% cytobands$chrom)
      msg <- c(msg, "seqnames not found in cytobands\n")
  }

  if (is.null(msg)) return(TRUE)
  message(msg)
  FALSE
}

.checkCoverage <- function(
  msg, coverage, linecol, type, annotation, annotcol, ylim
) {

  if (missing(coverage)) return(msg)

  ## Add checks for coverage. This can be a list of BigWigFileLists, or a
  ## single BigWigFileList
  if (is(coverage, "list")) {
    chk <- vapply(coverage, is, logical(1), class2 = "BigWigFileList")
    if (!all(chk))
      msg <- c(msg, "All elements of 'coverage' must be a BigWigFileList\n")
    if ("" %in% names(coverage))
      msg <- c(msg, "'coverage' must a be a named list\n")
  } else {
    if (!is(coverage, "BigWigFileList"))
      msg <- c(msg, "'coverage' should be a BigWigFileList or a list\n")
  }
  if (!is.null(linecol) & type == "l") {
    ## If we have BigWigFileList, we need a vector or list of colours
    ## the same length. If named, names must match
    if (is(coverage, "BigWigFileList")) {
      linecol <- unlist(linecol)
      if (length(linecol) != length(coverage))
        msg <- c(
          msg, "linecol must be the same length os the coverage tracks\n"
        )
      if (!is.null(names(linecol)) & !all(names(coverage) %in% names(linecol)))
        msg <- c(
          msg, "All elements of coverage must be named in linecol\n"
        )
    }
    ## If passing a list of BigWigFileList objects colours must be a list of
    ## the same length. Each element must also be NULL, or match
    if (is(coverage, "list")) {
      nm <- names(coverage)
      if (!is.list(linecol)) {
        msg <- c(msg, "linecol must be a named list\n")
      } else {
        if (!all(nm %in% names(linecol))) {
          msg <- c(msg, "All elements of coverage must be named in linecol\n")
        } else {
          cov_len <- vapply(coverage, length, integer(1))
          col_len <- vapply(linecol, length, integer(1))
          cov_len <- cov_len[names(which(col_len > 0))]
          if (!all(cov_len == col_len[names(cov_len)]))
            msg <- c(
              msg,
              "All elements of linecol must be NULL, or match coverage\n"
            )
        }
      }
    }
  }

  if (!is.null(ylim)) {
    if (is(coverage, "list")) {
      # If coverage is a list, ylim can be a vector passed to all elements
      # or a named list of numerics
      if (!is(ylim, "list")) {
        if (length(ylim) < 2 | !is.numeric(ylim))
          msg <- c(
            msg, "ylim can be a named list or numeric vector of length >= 2\n"
          )
      } else {
        if (!all(names(coverage) %in% names(ylim)))
          msg <- c(
            msg, "All elements of coverage must also be named in ylim\n"
          )
        if (!all(vapply(ylim, length, integer(1)) >= 2))
          msg <- c(
            msg, "All supplied elements of ylim must be of length >= 2\n"
          )
      }
    }
    if (is(coverage, "BigWigFileList")) {
      if (length(ylim) < 2 | !is.numeric(ylim))
        msg <- c(
          msg, "ylim should be passed as a numeric vector of  length >= 2"
        )
    }
  }

  if (missing(annotation)) return(msg)

  all_annot <- c()
  ## We need to check that it's a GRangesList if coverage is a BigWigFileList
  if (is(coverage, "BigWigFileList")) {
    if (!is(annotation, "GRangesList")) {
      msg <- c(msg, "annotation must be a GRangesList\n")
    } else {
      all_annot <- names(annotation)
    }
  }

  ## Or a list of GRangesLists if coverage is a list of BigWigFileLists
  if (is(coverage, "list")) {
    if (!is(annotation, "list")) {
      msg <- c(msg, "annotation must be a list of GRangesList objects\n")
    } else {
      is_grl <- vapply(annotation, is, logical(1), class2 = "GRangesList")
      if (!all(is_grl)) {
        msg <- c(msg, "annotation must be a list of GRangesList objects\n")
      } else {
        all_annot <- unlist(lapply(annotation, names))
      }
    }
  }

  if (missing(annotcol)) return(msg)

  ## Also we need to check that colours match if specified
  if (!missing(annotcol)) {
    all_cols <- names(annotcol)
    miss <- setdiff(all_annot, all_cols)
    if (length(miss) > 0)
      msg <- c(msg, "Colours not specified for ", miss, "\n")
  }

  msg

}

