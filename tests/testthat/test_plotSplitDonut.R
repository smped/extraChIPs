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
    expect_error(plotSplitDonut(df), 'argument .+ is missing')
    expect_error(plotSplitDonut(df, inner = "TF1"), 'argument "outer" is missing')
    expect_error(
        plotSplitDonut(df, inner = "TF1", outer = "TF1"),
        'inner != outer is not TRUE'
    )
    expect_error(
        plotSplitDonut(df, inner = "TF1", outer = "feature", scale_by = "a"),
        "all\\(c\\(inner, outer, scale_by\\) %in% colnames\\(object\\)\\) is not TRUE"
    )
    expect_error(
        plotSplitDonut(df, inner = "TF1", outer = "feature", scale_by = "TF1"),
        "is.numeric.+"
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
    expect_true(all(p[[1]]$data$angle == 0))
    expect_true(is(p[[2]], "wrapped_patch"))
    expect_true(is(p[[3]], "spacer"))
    expect_true(
        sum(dplyr::filter(p[[1]]$data, ring == "outer")$p) == 1
    )

    p <- plotSplitDonut(
        gr, inner = "TF1", outer = "feature", scale_by = "width",
        outer_glue = "{.data[[outer]]}\n({round(n, 1)}kb)",
        explode_inner = "Up", explode_outer = "Prom", explode_r = 0.2,
        inner_palette = hcl.colors(3, "Spectral", rev = TRUE),
        outer_palette = hcl.colors(3, "Cividis"),
        inner_rotate = TRUE, outer_rotate = TRUE
    )
    expect_true(is(p[[3]], "wrapped_patch"))
    expect_equal(sum(p[[1]]$data$n), sum(width(gr)/1e3)*2)
    expect_true(!any(p[[1]]$data$angle == 0))
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
    expect_equal(sum(grepl("kb", p[[1]]$data$lab)), 9)

    p <- plotSplitDonut(df, inner = "TF1", outer = "TF2", outer_p_by = "inner")
    expect_true(
        sum(dplyr::filter(p[[1]]$data, ring == "outer")$p) == 3
    )

})


test_that("DataFrame objects work as expected", {
    p <- plotSplitDonut(DataFrame(df), inner = "TF1", outer = "TF2")
    expect_true(is(p, "patchwork"))
    expect_true(is(p[[1]], "gg"))
    expect_true(is(p[[2]], "wrapped_patch"))
    expect_true(is(p[[3]], "spacer"))
})

test_that("label types are handled correctly",{
    p <- plotSplitDonut(df, inner = "TF1", outer = "TF2")
    p_build <- ggplot_build(p[[1]])
    expect_true(length(p_build$data) == 4)
    expect_true(
        all(vapply(p_build$data, is, logical(1), "data.frame"))
    )
    p <- plotSplitDonut(
        df, inner = "TF1", outer = "TF2", total_label = "none",
        inner_label = "none", outer_label = "none"
    )
    p_build <- ggplot_build(p[[1]])
    expect_true(length(p_build$data) == 1)
})
