gr <- GRanges(rep(c("chr1:1-10"), 2))
gr$id <- paste0("range", seq_along(gr))
gr$gene <- "gene1"

test_that("Default settings of distinct are consistent", {
  expect_equal(unique(granges(gr)), distinctMC(gr))
  d_gr <- distinctMC(gr, gene)
  expect_equal(d_gr$gene, "gene1")
  expect_equal(distinctMC(GRanges()), GRanges())
})

