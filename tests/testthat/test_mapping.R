# We need a combination of ranges which do & don't map
genes <- GRanges(c("chr1:1-10:*", "chr1:25-30:-", "chr1:31-40:+"))
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
gr$via_prom <- CharacterList(
  list("gene1", NULL, NULL, c("gene2", "gene3"), NULL, "gene1")
  )
gr$via_enh_gap0 <- CharacterList(
  list(NULL, NULL, NULL, NULL, NULL, NULL)
)
gr$via_enh_gap5 <- CharacterList(
  list(NULL, NULL, "gene2", NULL, NULL, NULL)
)
gr$via_gi <- CharacterList(
  list("gene1", NULL, NULL, NULL, NULL, "gene1")
)

test_that(".mapFeatures produces expected output from promoters", {
  mapped_prom <- .mapFeatures(gr, prom, genes, "gene_id", 0, 0)
  expect_true(is(mapped_prom, "tbl_df"))
  expect_equal(
    table(mapped_prom$range),
    structure(1:2, .Dim = 2L, .Dimnames = structure(list(c("chr1:0", "chr1:30")), .Names = ""), class = "table")
  )
  expect_equal(mapped_prom$gene_id, c("gene1", "gene2", "gene3"))
})

test_that(".mapFeatures produces expected output from enhancers", {
  expect_equal(
    .mapFeatures(gr, enh, genes, "gene_id", 0, 0)$gene_id, rep(NA_character_, 2)
  )
  expect_equal(
    .mapFeatures(gr, enh, genes, "gene_id", 0, 5)$gene_id,
    c("gene2", NA)
  )
})

test_that(".mapFeatures returns NULL & errors where expected", {
  expect_null(.mapFeatures(NULL))
  expect_null(.mapFeatures(gr, NULL))
  expect_null(.mapFeatures(gr))
  expect_error(.mapFeatures(gr, prom, genes, "", 0, 0))
})

# test_that(".mapGi produces the expected output",{
#
# })

test_that(".mapWithin returns NULL & errors when expected", {
  expect_null(.mapWithin(NULL))
  expect_null(.mapWithin(gr))
  expect_error(.mapWithin(gr, genes, ""))
})
