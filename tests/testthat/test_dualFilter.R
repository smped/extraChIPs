## These are extremeley rudimentary tests, just probing the input classes
## More detailed tests checking the output consistency are still required
## Detailed tests for bam files are also required

test_that("Erroneous input specs are caught", {

    rse <- SummarizedExperiment(rowRanges = GRanges())
    rse2 <- SummarizedExperiment(rowRanges = GRanges("chr1:1"))
    expect_error(dualFilter(""))
    expect_error(dualFilter(rse, ""))
    expect_error(dualFilter(rse, rse, ""))
    expect_error(dualFilter(rse, rse, GRanges(), ""))
    expect_error(dualFilter(rse, rse, GRanges(), 0.5, ""))
    expect_error(dualFilter(rse, rse2, GRanges(), 0.5, TRUE))
    expect_error(dualFilter(rse, rse, GRanges(), 0.5, TRUE))

})
