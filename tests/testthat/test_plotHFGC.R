gr1 <- GRanges("chr1:1-10")
gr2 <- GRanges("chr1:21-30")
hic <- InteractionSet::GInteractions(gr1, gr2)
hiccol <- list(anchors = "lightblue", interactions = "red")
feat <- GRangesList(a = GRanges("chr1:8"))
genes <- gr1
mcols(genes) <- DataFrame(
  feature = "exon", gene = "ENSG", exon = "ENSE", transcript = "ENST",
  symbol = "ID"
)
cyto_df <- data.frame(
  chrom = "chr1", chromStart = 1, chromEnd = 100, name = "p1", gieStain = "gneg"
)

test_path <- system.file("tests", "test.bw", package = "rtracklayer")
test_bw <- BigWigFileList(test_path)
names(test_bw) <- "a"

test_that(".checkHFGCArgs catches GRanges issues", {
  expect_error(
    .checkHFGCArgs(zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l")
  )
  expect_message(
    .checkHFGCArgs(
      gr = "", zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l"
    ),
    "must be a GRanges object"
  )
  expect_message(
    .checkHFGCArgs(
      gr = GRanges(), zoom = 1, shift = 0, max = 1e7, axistrack = TRUE,
      type = "l"
    ),
    "Cannot be an empty range"
  )
  expect_true(
    .checkHFGCArgs(
      gr = gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l"
    )
  )
})

test_that("Malformed numerics/bools are caught", {
  expect_message(
    .checkHFGCArgs(
      gr1, zoom = "", shift = 0, max = 1e7, axistrack = TRUE, type = "l"
    ),
    "zoom/shift/max must be numeric or coercible to numeric"
  )
  expect_message(
    .checkHFGCArgs(
      gr1, zoom = 1, shift = 0, max = 1e7, axistrack = "", type = "l"
    ),
    "axistrack must be logical"
  )
})

test_that("Malformed HiC Interactions are caught", {
  expect_true(
    .checkHFGCArgs(
      gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
      hic = hic, hiccol = hiccol
    )
  )
  expect_message(
    .checkHFGCArgs(
      gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
      hic = c(), hiccol = hiccol
    ),
    "'hic' must be provided as a GInteractions object"
  )
  expect_message(
    .checkHFGCArgs(
      gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
      hic = hic, hiccol = hiccol[1]
    ),
    "'hiccols' must be a list with name: interactions"
  )

})

test_that("Malformed features are caught", {
  expect_true(
    .checkHFGCArgs(
      gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
      features = GRangesList()
    )
  )
  expect_message(
    .checkHFGCArgs(
      gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
      features = c()
    ),
    "'features' must be provided as a GRangesList"
  )
  expect_message(
    .checkHFGCArgs(
      gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
      features = GRangesList(a = GRanges(), GRanges())
    ),
    "All elements of 'features' must be explicitly named"
  )
  expect_message(
    .checkHFGCArgs(
      gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
      features = GRangesList(a = GRanges()), featurecol = c(b = "blue")
    ),
    "All elements of 'features' must be in 'featurecol'"
  )
})

test_that("Malformed genes are caught", {
  expect_true(
    .checkHFGCArgs(
      gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
      genes = genes
    )
  )
  expect_message(
    .checkHFGCArgs(
      gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
      genes = GRanges()
    ),
    "'genes' must have an mcols component with the columns"
  )
  expect_message(
    .checkHFGCArgs(
      gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
      genes = GRangesList(a = gr1), genecol = "blue"
    ),
    "'genes' must have an mcols component with the columns"
  )
})

test_that("Malformed coverage parameters are caught", {
  expect_true(
    .checkHFGCArgs(
      gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
      coverage = test_bw, linecol = c(a = "blue")
    )
  )
  expect_true(
    .checkHFGCArgs(
      gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
      coverage = list(a = test_bw), linecol = c(a = "blue")
    )
  )
  expect_message(
    .checkHFGCArgs(
      gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
      coverage = list(a = c())
    ),
   "All elements of 'coverage' must be a BigWigFileList"
  )
  expect_message(
    .checkHFGCArgs(
      gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
      coverage = c()
    ),
    "'coverage' should be a BigWigFileList or a list"
  )
  expect_message(
    .checkHFGCArgs(
      gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
      coverage = list(a = test_bw), linecol =  "blue"
    ),
    "'linecol' must be a named vector of colours"
  )

})

test_that("Malformed cytobands are caught", {
  expect_true(
    .checkHFGCArgs(
      gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
      cytobands = cyto_df
    )
  )
  expect_message(
    .checkHFGCArgs(
      gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
      cytobands = cyto_df[1]
    ),
    "Columns in cytobands must be exactly"
  )
  expect_message(
    .checkHFGCArgs(
      gr1, zoom = 1, shift = 0, max = 1e7, axistrack = TRUE, type = "l",
      cytobands = cyto_df[c(),]
    ),
    "seqnames not found in cytobands"
  )
})

test_that("Ideogram forms correctly", {
  ideo <- .makeIdeoTrack(gr1, cyto_df, 12)
  expect_true(is(ideo, "IdeogramTrack"))
  expect_message(
    .makeIdeoTrack(gr1, .fontsize = 12),
    "Could not find cytogenetic bands for genome NA"
  )
})

test_that("HiC Track forms correctly", {
  hic_track <- .makeHiCTrack(
    hic, gr1, .tracksize = 2, .fontsize = 10, .cex = 1, .rot = 90, .col = hiccol
  )
  expect_true(is(hic_track, "InteractionTrack"))
  expect_equal(hic_track@dp@pars$col.anchors.line, hiccol$anchors)
  expect_equal(hic_track@dp@pars$col.anchors.fill, hiccol$anchors)
  expect_equal(hic_track@dp@pars$col.interactions, hiccol$interactions)
  expect_equal(hic_track@dp@pars$fontsize, 10)
  expect_equal(hic_track@dp@pars$size, 2)
  expect_equal(hic_track@dp@pars$rotation.title, 90)
  expect_null(.makeHiCTrack(hic, GRanges()))
})

test_that("Feature Track forms correctly", {
  feat_track <- .makeFeatureTrack(
    feat, gr1, .fontsize = 12, list(a = "blue"), .cex = 1, .tracksize = 1,
    .rot = 0
  )
  expect_true(is(feat_track, "AnnotationTrack"))
  expect_equal(feat_track@dp@pars$fill, c(a = "blue"))
  expect_null(
    .makeFeatureTrack(
      feat, gr2, .fontsize = 12, list(a = "blue"), .cex = 1, .tracksize = 1,
      .rot = 0
    )
  )
})

test_that("Genes Track forms correctly", {
  gene_track <- .makeGeneTracks(
    genes, gr1, "meta", .tracksize = 1, .cex = 1, .rot = 0, .fontsize = 12
  )
  expect_true(is(gene_track, "GeneRegionTrack"))
  expect_equal(length(gene_track), 1)
  expect_equal(gene_track@dp@pars$fill, "#FFD58A")
  expect_equal(gene_track@dp@pars$fontsize, 12)
  expect_equal(gene_track@dp@pars$cex, 1)
  expect_equal(gene_track@dp@pars$collapseTranscripts, "meta")
  expect_equal(gene_track@dp@pars$rotation, 0)
  expect_equal(
    .makeGeneTracks(genes, gr2, "meta", .tracksize = 1, .cex = 1, .rot = 0),
    list(NULL)
  )

  gene_track <- .makeGeneTracks(
    GRangesList(a = genes),
    gr1, "meta", .tracksize = 1, .cex = 1, .rot = 0, .fontsize = 12
  )
  expect_true(is(gene_track, "list"))
  expect_true(is(gene_track[[1]], "GeneRegionTrack"))
  expect_equal(names(gene_track[[1]]), "A")
  expect_equal(gene_track[[1]]@dp@pars$fill, "#E41A1C")

})

test_that("Coverage Track forms correctly", {
  cov_track <- .makeCoverageTracks(
    .coverage = test_bw, .gr = GRanges("chr2:1-10"),
    .fontsize = 12, .type = "l", .gradient = "blue", .tracksize = 1, .cex = 1,
    .rot = 0, .linecol = c(a = "red")
  )
  expect_true(is(cov_track, "list"))
  expect_true(is(cov_track[[1]], "DataTrack"))
  expect_equal(length(cov_track), 1)
  expect_equal(cov_track[[1]]@dp@pars$gradient, "blue")
  expect_equal(cov_track[[1]]@dp@pars$size, 1)
  expect_equal(cov_track[[1]]@dp@pars$rotation, 0)
  expect_equal(cov_track[[1]]@dp@pars$col, c(a = "red"))
  expect_equal(levels(cov_track[[1]]@dp@pars$groups), "a")
  expect_equal(
    suppressWarnings(.makeCoverageTracks(.coverage = test_bw, .gr = gr1)),
    list(NULL)
  )
  expect_warning(
    .makeCoverageTracks(.coverage = test_bw, .gr = gr1),
    "'which' contains seqnames not known to BigWig file: chr1"
  )

})
