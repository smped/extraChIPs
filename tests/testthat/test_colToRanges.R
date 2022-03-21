sq <- Seqinfo("chr1", 100, FALSE, "test")
x <- GRanges(c("chr1:1-10", "chr1:6-15", "chr1:51-60"), seqinfo = sq)
df <- data.frame(logFC = rnorm(3), logCPM = rnorm(3,8), p = 10^-rexp(3))
gr <-mergeByCol(x, df, col = "logCPM", pval = "p")

test_that("Coercion is correct where it should be", {
  new_gr <- colToRanges(gr, "keyval_range")
  expect_s4_class(new_gr, "GRanges")
  expect_equal(
    colnames(mcols(new_gr)),
    c("n_windows", "n_up", "n_down", "logCPM", "logFC", "p", "p_fdr")
  )
  expect_equal(length(new_gr), 2)
  df$gr <- as.character(x)
  new_gr <- colToRanges(df, "gr", seqinfo = sq)
  expect_s4_class(new_gr, "GRanges")
  expect_equal(seqinfo(new_gr), sq)
})

test_that("Coercion fails where it should", {
  expect_error(colToRanges(x, ""))
  expect_error(colToRanges(df, "logFC"))
})
