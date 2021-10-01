test_that("drop.existing behaves correctly", {
    gr <- GRanges("chr1:1-10")
    gr$p <- 0.2
    df <- as_tibble(gr)
    expect_true(is(df, "tbl_df"))
    expect_equal(c("range", "p"), colnames(df))
    expect_equal(ncol(as_tibble(gr, rangeAsChar = FALSE)), 6)
})

