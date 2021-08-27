test_that("drop.existing behaves correctly", {
    gr <- GRanges("chr1:1-10")
    gr$p <- 0.2
    df <- rangesAsTibble(gr)
    expect_error(rangesAsTibble(gr, "p"))
    expect_error(rangesAsTibble(df))
    expect_true(is(df, "tbl_df"))
    expect_true(all(c("range", "p") %in% colnames(df)))
})

test_that("mcolsAsTibble behaves correctly", {
    gr <- GRanges("chr1:1-10")
    gr$p <- 0.2
    df <- mcolsAsTibble(gr)
    expect_true(colnames(df) == "p")
    expect_error(mcolsAsTibble(df))
})
