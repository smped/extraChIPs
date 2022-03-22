gr <- GRanges(rep(c("chr1:1-10"), 2))
gr$id <- paste0("range", seq_along(gr))
gr$gene <- "gene1"

test_that("Default settings for chopMC work", {
  chop_gr <- chopMC(gr)
  expect_equal(reduce(gr), granges(chop_gr))
  expect_equal(chop_gr$gene, "gene1")
  expect_equal(chop_gr$id[[1]], gr$id)
  chop_gr <- chopMC(gr, simplify = FALSE)
  expect_true(is(chop_gr$gene, "CharacterList"))
})
