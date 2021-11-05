# We need a combination of ranges which do & don't map
genes <- GRanges(c("chr1:2-10:*", "chr1:25-30:-", "chr1:31-40:+"))
genes$gene_id <- paste0("gene", seq_along(genes))
## Add a promoter for each gene
prom <- promoters(genes, upstream = 1, downstream = 1)
## Add an enhancer in gene2 & to the right of gene3
enh <- GRanges(c("chr1:20", "chr1:50"))
enh$which <- c("gene2", "gene1 (gi)")
## Map the last enhancer to gene1's promoter
gi <- InteractionSet::GInteractions(prom[1], enh[2])
mcols(gi) <- c()
## Now define key ranges
gr <- GRanges(paste0("chr1:", seq(0, 50, by = 10)))
gr$via_ol <- CharacterList(
  list(NULL, "gene1", NULL, c("gene2", "gene3"), "gene3", NULL)
)
gr$via_prom <- CharacterList(
  list("gene1", NULL, NULL, c("gene2", "gene3"), NULL, NULL)
  )
gr$via_enh_gap5 <- CharacterList(
  list(NULL, NULL, "gene2", NULL, NULL, NULL)
)
gr$via_gi <- CharacterList(
  list("gene1", NULL, NULL, NULL, NULL, "gene1")
)

test_that(".mapFeatures produces expected output from promoters", {
  mapped_prom <- .mapFeatures(gr, prom, genes, "gene_id", 0, 0)
  expect_true(is(mapped_prom, "data.frame"))
  expect_equal(
    table(mapped_prom$range),
    structure(
      1:2, .Dim = 2L,
      .Dimnames = structure(list(c("chr1:0", "chr1:30")), .Names = ""),
      class = "table"
    )
  )
  expect_equal(mapped_prom$gene_id, c("gene1", "gene2", "gene3"))
})

test_that(".mapFeatures produces expected output from enhancers", {
  expect_equal(.mapFeatures(gr, enh, genes, "gene_id", 0, 5)$gene_id, "gene2")
})

test_that(".mapFeatures returns NULL & errors where expected", {
  expect_null(.mapFeatures(gr, enh, genes, "gene_id", 0, 0))
  expect_null(.mapFeatures(NULL))
  expect_null(.mapFeatures(gr, NULL))
  expect_null(.mapFeatures(gr))
  expect_error(.mapFeatures(gr, prom, genes, "", 0, 0))
})

test_that(".mapGi produces the expected output",{
  mapped_gi <- .mapGi(gr, gi, genes, "gene_id", 0, 0)
  expect_true(is(mapped_gi, "data.frame"))
  expect_equal(
    table(mapped_gi$range),
    structure(
      c(`chr1:0` = 1L, `chr1:50` = 1L), .Dim = 2L,
      .Dimnames = structure(list(c("chr1:0", "chr1:50")), .Names = ""),
      class = "table"
    )
  )
  expect_equal(mapped_gi$gene_id, c("gene1", "gene1"))
})

test_that("mapGi returns NULL & errors where expected", {
  expect_error(.mapGi(gr, gi, genes, "", 0, 0))
  expect_null(.mapGi(gr[NULL], gi, genes, "gene_id", 0, 0))
  expect_null(.mapGi(gr, gi[NULL], genes, "gene_id", 0, 0))
})

test_that(".mapWithin produces the expected output", {
  mapped_within <- .mapWithin(gr, genes, "gene_id", 0)
  expect_true(is(mapped_within, "data.frame"))
  expect_equal(
    table(mapped_within$range),
    structure(
      c(`chr1:10` = 1L, `chr1:30` = 2L, `chr1:40` = 1L), .Dim = 3L,
      .Dimnames = structure(list(c("chr1:10", "chr1:30", "chr1:40")), .Names = ""),
      class = "table"
    )
  )
  expect_equal(mapped_within$gene_id, paste0("gene", c(1, 2, 3, 3)))
})

test_that(".mapWithin returns NULL & errors when expected", {
  expect_null(.mapWithin(NULL))
  expect_null(.mapWithin(gr))
  expect_error(.mapWithin(gr, genes, ""))
})

test_that("mapByFeature produces the correct output", {
  ol_only <- mapByFeature(gr, genes, cols = "gene_id", gr2gene = 0)$gene_id
  expected <- new(
    "CompressedCharacterList", elementType = "character",
    elementMetadata = NULL, metadata = list(),
    unlistData = c("gene1", "gene2", "gene3", "gene3"),
    partitioning = new(
      "PartitioningByEnd", end = c(0L, 1L, 1L, 3L, 4L, 4L),
      NAMES = NULL, elementType = "ANY", elementMetadata = NULL, metadata = list()
    )
  )
  expect_equal(ol_only, expected)
  ## Update promoter only mappings
  expected[[1]] <- "gene1"
  expect_equal(
    mapByFeature(gr, genes, prom = prom, cols = "gene_id", gr2gene = 0)$gene_id,
    expected
  )
  ## Update enhancer only mappings
  expected[[3]] <- "gene2"
  expect_equal(
    mapByFeature(
      gr, genes, prom = prom, enh = enh, cols = "gene_id",
      gr2gene = 0, enh2gene = 5
    )$gene_id,
    expected
  )
  ## Update interaction mappings
  expected[[6]] <- "gene1"
  expect_equal(
    mapByFeature(
      gr, genes, prom = prom, enh = enh, gi = gi,
      cols = "gene_id", gr2gene = 0, enh2gene = 5
    )$gene_id,
    expected
  )

})

test_that("mapByFeature errors correctly", {
  expect_error(mapByFeature(NULL))
  expect_error(suppressMessages(mapByFeature(gr, genes, cols = "")))
})

test_that("mappingArg Checks behave correctly", {
  expect_error(.checkMappingArgs(), "query ranges must be provided")
  expect_error(
    .checkMappingArgs(NULL), "query object must be a GenomicRanges object"
  )
  expect_error(
    .checkMappingArgs(GRanges(), .cols = NULL),
    "At least one gene/ID column must be specified"
  )
  expect_message(
    .checkMappingArgs(gr, genes, prom, enh, gi, "", 0),
    "No requested columns were able to be found"
  )
  expect_message(
    .checkMappingArgs(gr, genes, prom, enh, gi, "gene_id", ""),
    "All distance arguments must be numeric"
  )
})