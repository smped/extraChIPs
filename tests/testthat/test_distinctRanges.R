library(tidyselect)
gr <- GRanges(rep(c("chr1:1-10"), 2))
gr$id <- paste0("range", seq_along(gr))
gr$gene <- "gene1"

test_that("Default settings of distinctRanges are consistent", {
  expect_equal(gr, distinctRanges(gr))
  d_gr <- distinctRanges(gr, all_of("gene"))
  expect_equal(d_gr$gene, "gene1")
})

