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

test_that("mapByFeature produces the correct output", {
  ol_only <- mapByFeature(gr, genes, cols = "gene_id", gr2gene = 0)$gene_id
  expected <- new(
    "CompressedCharacterList", elementType = "character",
    elementMetadata = NULL, metadata = list(),
    unlistData = c("gene1", "gene2", "gene3"),
    partitioning = new(
      "PartitioningByEnd", end = c(0L, 1L, 1L, 2L, 3L, 3L),
      NAMES = NULL, elementType = "ANY", elementMetadata = NULL, metadata = list()
    )
  )
  expect_equal(ol_only, expected)
  ## Update promoter only mappings
  expected[[1]] <- "gene1"; expected[[4]] <- c("gene2", "gene3")
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

