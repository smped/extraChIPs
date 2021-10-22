gr <- GRanges("chr1:1-10")
prom <- GRanges("chr1:1")
enh <- GRanges("chr1:10")
hic <- InteractionSet::GInteractions(gr, GRanges("chr1:51-60"))
hic$gene_id <- list(c("a", "b"))
genes <- GRanges(c("chr1:5", "chr1:55"))
genes$gene_id <- c("a", "b")

test_that("checkArgs catches everything", {
  expect_true(.checkMappingArgs(gr, prom, enh, hic, genes, "gene_id"))
  expect_message(
    .checkMappingArgs(gr, prom, enh, hic, genes, NULL),
    "No columns found in 'cols' for obtaining gene information"
  )
  expect_message(
    .checkMappingArgs(NULL, NULL, NULL, hic, NULL, "gene_id"),
    "'genes'.+'gr'.+'prom'.+'enh' must be a GRanges object",
  )
  expect_message(
    .checkMappingArgs(gr, prom, enh, NULL, genes, "gene_id"),
    "'hic' must be provided as GInteractions"
  )
  expect_message(
    .checkMappingArgs(gr, prom, enh, hic, granges(genes), "gene_id"),
    "All valid columns in 'cols' must be in both 'hic' and 'genes'"
    )
  expect_error(.checkMappingArgs(gr, prom, enh, .cols = "gene_id"))
  expect_error(
    suppressMessages(rangesToGenes(NULL, prom, enh, hic, genes, "gene_id"))
  )
})
