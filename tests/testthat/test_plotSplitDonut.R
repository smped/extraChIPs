set.seed(200)
df <- data.frame(
    feature = sample(c("Promoter", "Enhancer", "Intergenic"), 200, replace = TRUE),
    TF1 = sample(c("Up", "Down", "Unchanged"), 200, replace = TRUE),
    TF2 = sample(c("Up", "Down", "Unchanged"), 200, replace = TRUE)
)
gr <- makeGRangesFromDataFrame(
    tibble(
        chr = "chr1", start = sample.int(1e6, 200), end = start + rpois(200, 100)
    )
)
mcols(gr) <- df

test_that("plotSplitDonut Errors where expected", {

    expect_error(plotSplitDonut(NULL))
    expect_error(plotSplitDonut(df), 'argument "inner" is missing')
    expect_error(plotSplitDonut(df, inner = "TF1"), 'argument "outer" is missing')
    expect_error(
        plotSplitDonut(df, inner = "TF1", outer = "TF1"),
        'inner != outer is not TRUE'
    )
    expect_error(
        plotSplitDonut(df, inner = "TF1", outer = "feature", scale_by = ""),
        "all\\(c\\(inner, outer, scale_by\\) %in% colnames\\(object\\)\\) is not TRUE"
    )
    expect_error(
        plotSplitDonut(df, inner = "TF1", outer = "feature", scale_by = "TF1"),
        "is.numeric\\(object\\[\\[scale_by\\]\\]\\) is not TRUE"
    )

})

test_that("plotSplitDonut produces the expected output", {

    p <- plotSplitDonut(df, inner = "TF1", outer = "TF2")
    expect_equal(dim(p[[1]]$data), c(12, 19))
    expect_equal(sum(p[[1]]$data$n), 200 * 2)
    expect_equal(
        levels(p[[1]]$data$colour), paste("TF1", c("Down", "Unchanged", "Up"))
    )
    expect_true(all(!p[[1]]$data$explode))
    expect_true(is(p[[2]], "wrapped_patch"))
    expect_true(is(p[[3]], "spacer"))

    p <- plotSplitDonut(
        gr, inner = "TF1", outer = "feature", scale_by = "width",
        explode_inner = "Up", explode_outer = "Prom", explode_r = 0.2,
        inner_palette = hcl.colors(3, "Spectral", rev = TRUE),
        outer_palette = hcl.colors(3, "Cividis"), cat_outer = FALSE
    )
    expect_true(is(p[[3]], "wrapped_patch"))
    expect_equal(sum(p[[1]]$data$n), sum(width(gr)/1e3)*2)
    expect_equal(
        levels(p[[1]]$data$colour),
        c(
            paste("TF1", sort(unique(gr$TF1))),
            paste("feature", sort(unique(gr$feature)))
        )
    )
    expect_equal(which(p[[1]]$data$explode), 12)
    expect_equal(p[[1]]$data$x[[12]], 1.7)
    expect_true(all(p[[1]]$data$alpha == 1))
    expect_true(all(!grepl("feature", p[[1]]$data$lab)))

})


test_that("DataFrame objects work as expected", {
    p <- plotSplitDonut(DataFrame(df), inner = "TF1", outer = "TF2")
    expect_true(is(p, "patchwork"))
    expect_true(is(p[[1]], "gg"))
    expect_true(is(p[[2]], "wrapped_patch"))
    expect_true(is(p[[3]], "spacer"))
})
