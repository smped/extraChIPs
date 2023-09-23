data('peaks')
names(peaks) <- gsub("_peaks.+", "", names(peaks))
df <- data.frame(sample = names(peaks), treat = rep(c("A", "B"), each = 3))

test_that("General plots work", {
    p <- plotGrlCol(peaks)
    expect_true(is(p, "gg"))
    p <- plotGrlCol(peaks, var = "signal")
    expect_true(is(p, "gg"))
    p <- plotGrlCol(peaks, df = df, fill = treat)
    expect_true(is(p, "gg"))
})

test_that("Errors are caught", {
    expect_error(plotGrlCol(peaks, var = ""))
    expect_error(plotGrlCol(peaks, df = df, .id = ""))
    expect_error(plotGrlCol(peaks, df = df, fill = ""))
    expect_error(plotGrlCol(peaks, df = df, colour = ""))
})
