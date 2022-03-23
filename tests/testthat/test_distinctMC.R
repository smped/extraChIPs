gr <- GRanges(rep(c("chr1:1-10"), 2))
gr$id <- paste0("range", seq_along(gr))
gr$gene <- "gene1"

test_that("distinctMC errors when expected", {
  expect_error(distinctMC(gr, ""), "Requested columns absent")
})

test_that("Default settings of distinctMC are consistent", {
  expect_equal(gr, distinctMC(gr))
  d_gr <- distinctMC(gr, "gene")
  expect_equal(d_gr$gene, "gene1")
  expect_equal(distinctMC(granges(gr)), unique(granges(gr)))
  expect_equal(gr, distinctMC(gr)) # Notthing should happen for the test data
})

