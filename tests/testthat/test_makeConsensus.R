a <- GRanges("chr1:11-20")
a$score <- 1
b <- GRanges(c("chr1:18-22", "chr1:1-5"))
b$score <- c(0.6, 0.3)
grl <- GRangesList(a = a, b = b)

test_that("makeConsensus errors correctly", {
    expect_error(makeConsensus(as.list(grl)))
    expect_error(makeConsensus(grl, var = ""), "Couldn't find column")
})

test_that("makeConsensus returns correct output", {
    gr <- makeConsensus(grl)
    expect_length(gr, 2)
    expect_equal(colnames(mcols(gr)), c("a", "b", "n"))
    expect_length(makeConsensus(grl, p = 1), 1)
    gr <- makeConsensus(grl, var = "score")
    expect_true(is(gr$score, "CompressedNumericList"))
})
