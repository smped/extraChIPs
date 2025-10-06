ex <- list(
    x = letters[1:5], y = letters[c(6:15, 26)], z = letters[c(2, 10:25)]
)
set.seed(100)
grl <- GRangesList(
    a = GRanges(c("chr1:1-10", "chr1:21-30", "chr1:31-40")),
    b = GRanges(c("chr1:12-15", "chr1:21-30", "chr1:46-50"))
)
grl$a$score <- rnorm(3)
grl$b$score <- rnorm(3)

test_that(".mcolnames works",{
    expect_equal(.mcolnames(grl[[1]]), "score")
})

test_that("plotSingleVenn produces correct output", {
    p <- .plotSingleVenn(ex[1])
    expect_equal(length(p), 4)
    expect_equal(p[[1]]$params$x, 0.5)
    expect_equal(p[[1]]$params$y, 0.5)
    expect_equal(p[[3]]$label, "5")
    expect_error(.plotSingleVenn(ex))
})

test_that("plotDoubleVenn produces correct output", {
    p <- .plotDoubleVenn(ex[1:2])
    expect_equal(length(p), 8)
    expect_equal(p[[1]]$params$x, 0.34, tolerance = 1e-2)
    expect_equal(p[[2]]$params$x, 0.725, tolerance = 1e-2)
    expect_equal(p[[1]]$params$y, 0.5)
    expect_equal(p[[2]]$params$y, 0.5)
    expect_equal(p[[5]]$label, "11")
    expect_equal(p[[6]]$label, "5")
    expect_error(.plotDoubleVenn(ex))
})

test_that("plotTripleVenn produces correct output", {
    p <- .plotTripleVenn(ex)
    expect_equal(length(p), 14)
    expect_equal(
        vapply(p[1:3], function(x) x$params$x, numeric(1)), c(1:3)/4
    )
    expect_equal(
        vapply(p[1:3], function(x) x$params$y, numeric(1)), rep(0.5, 3)
    )
    expect_equal(
        vapply(p[7:14], function(x) x$label, character(1)),
        c(5, 6, 10, 1, 4, "y", "z", "x")
    )
    expect_error(.plotTribleVenn(ex[1:2]))
})

test_that("plotOverlaps dispatches type = 'auto' correctly", {
    p <- plotOverlaps(ex)
    expect_equal(length(p), 14)
    p <- plotOverlaps(ex, type = 'venn')
    expect_equal(length(p), 14)
    p <- plotOverlaps(ex, type = 'upset')
    expect_true(is(p, 'patchwork'))
    expect_equal(
        p@data,
        structure(
            list(
                intersect = structure(
                    c(4L, 5L, 5L, 3L, 2L, 2L, 1L), levels = c("1", "2", "3", "4", "5"),  class = "factor"
                ),
                degree = structure(
                    c(1L, 2L, 2L, 1L, 2L, 2L, 1L), levels = c("1", "2"), class = "factor"
                ),
                set = structure(
                    c(1L, 1L, 3L, 2L, 2L, 3L, 3L), levels = c("x", "y", "z"),  class = "factor"
                ),
                in_group = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
                set_int = c(1L, 1L, 3L, 2L, 2L, 3L, 3L),
                x = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
                y = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE),
                z = c(FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE)
            ),
            row.names = c(NA, -7L),
            class = c("tbl_df", "tbl", "data.frame")
        )
    )
    expect_equal(length(p@layers), 4)
    expect_error(
        plotOverlaps(ex[1], type = "upset"),
        "UpSet plots can only be drawn using more than one group"
    )
    # expect_equal(length(p$patches$plots), 3)
})

test_that("plotOverlaps adds annotations as expected", {
    p <- plotOverlaps(grl, type = 'upset', var = 'score', set_col = "red")
    expect_true(is(p, 'patchwork'))
    expect_equal(length(p), 4)
    expect_equal(length(p[[2]]), 2)
    bp <- p[[2]][[1]]
    expect_true(is_ggplot(bp))
    expect_equal(
        bp@data$range, c("chr1:1-10", "chr1:12-15", "chr1:21-40", "chr1:46-50")
    )
    expect_equal(
        bp@data$score,
        c(-0.502192350531457, 0.886784809417845, 0.0565284486730913, 0.318630087617032)
    )
    expect_equal(as_label(bp@mapping$y), "score")
    expect_true(is(bp@layers[[1]]$geom, "GeomBoxplot"))
    expect_equal(as_label(p[[3]]@layers$geom_bar$mapping$fill), "set")
    expect_true(is(p[[3]]@scales$scales[[3]], "ScaleDiscrete"))
    expect_equal(ggplot_build(p[[3]])@data[[2]]$fill, rep("red", 2))

})

test_that("GRL Input is handled as expect without var", {
    expect_true(is(plotOverlaps(grl), "gList"))
})

test_that("Simple errors are caught", {
    expect_error(
        plotOverlaps(grl, type = "upset", var = "x"), "'arg' should be.+"
    )
    grl2 <- GRangesList(
        lapply(grl, function(x) {x$letters <- c("a", "b", "c"); x})
    )
    expect_error(
        plotOverlaps(grl2, type = "upset", var = "letters"),
        "letters must contain numeric values"
    )
    expect_error(
        plotOverlaps(grl[1], type = "upset"), "UpSet plots can only.+"
    )
})
