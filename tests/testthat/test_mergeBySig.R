sq <- Seqinfo("chr1", 100, FALSE, "test")
x <- GRanges(c("chr1:1-10", "chr1:6-15", "chr1:51-60"), seqinfo = sq)
set.seed(1001)
df <- DataFrame(logFC = rnorm(3), logCPM = rnorm(3,8), p = rexp(3, 10))
mcols(x) <- df

test_that("combineTestc behaves correctly for GRanges",{
    new_gr <- mergeBySig(x, df, pval = "p")
    expect_equal(length(new_gr), 2)
    expect_equal(sum(new_gr$p %in% df$p), 2)
    expect_equal(seqinfo(new_gr$keyval_range), sq)
    expect_equal(length(mergeBySig(x, pval = "p")), 2)
    expect_s4_class(mergeBySig(x, df, pval = "p")$keyval_range, "GRanges")
    expect_equal(as.character(new_gr$keyval_range), as.character(x)[-2])
})

test_that("getBestTest behaves correctly",{
    new_gr <- mergeBySig(x, pval = "p", method = "best")
    expect_equal(length(new_gr), 2)
    expect_equal(sum(new_gr$p %in% df$p), 1)
    expect_equal(as.character(new_gr$keyval_range), as.character(x)[2:3])
})

test_that("minimalTest behaves correctly", {
    new_gr <- mergeBySig(x, pval = "p", method = "min")
    expect_equal(sum(new_gr$p %in% df$p), 1)
    expect_equal(as.character(new_gr$keyval_range), as.character(x)[-2])
})

test_that("Errors appear where expected", {
    expect_error(mergeBySig(x, df[-1,], pval = "p"))
    expect_error(mergeBySig(x, df, min.sig.n = 3))
})

test_that("Weights are passed correctly", {
    gr_w <- mergeBySig(x, pval = "p", weights = c(1, 10, 1))
    expect_equal(as.character(gr_w$keyval_range), as.character(x)[2:3])
})
