x <- GRanges(c("chr1:1-10", "chr1:6-15", "chr1:51-60"))
df <- data.frame(logFC = rnorm(3), logCPM = rnorm(3,8), p = 10^-rexp(3))
mergeByCol(x, df, col = "logCPM", pval = "p")

test_that("Function behaves correctly for GRanges",{
  expect_equal(
    length(mergeByCol(x, df, col = "logCPM", pval = "p")), 2
  )
  mcols(x) <- df
  expect_equal(
    length(mergeByCol(x, col = "logCPM", pval = "p")), 2
  )
})

test_that("Function errors as expected", {
  expect_error(mergeByCol(x, df[1,], col = "logCPM", pval = "p"))
  expect_error(mergeByCol(x, df))
  expect_error(mergeByCol(x, df, col = "logFC"))
  expect_error(mergeByCol(x, df, col = ""))
})

test_that("P-value adjustment columns are correct", {
  expect_true(
    "fdr" %in% colnames(
      mcols(mergeByCol(x, df, col = "logCPM", pval = "p"))
    )
  )
  expect_equal(
    ncol(
      mcols(mergeByCol(x, df, col = "logCPM", pval = "p", p_adj_method = "none"))
    ),
    6
  )
})
