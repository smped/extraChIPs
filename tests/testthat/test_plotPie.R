set.seed(200)
df <- data.frame(
    feature = sample(c("Promoter", "Enhancer", "Intergenic"), 200, replace = TRUE),
    TF1 = sample(c("Up", "Down", "Unchanged"), 200, replace = TRUE),
    TF2 = sample(c("Up", "Down", "Unchanged"), 200, replace = TRUE)
)
data("ex_prom")

test_that("plotPie Errors where expected", {

    expect_error(plotPie(NULL))
    expect_error(plotPie(df), "The initial category must be defined as")
    expect_error(plotPie(df, fill = ""))
    expect_error(plotPie(df, fill = "feature", x = ""))
    expect_error(plotPie(df, fill = "feature", x = "TF1", y = ""))
    expect_error(plotPie(df, fill = "feature", scale_by = ""))

})

test_that(".plotSinglePie creates expected data structures", {

    p <- plotPie(df, "feature")
    expect_true(is(p, "gg"))
    expect_true(is.factor(p$data$feature))
    expect_equal(length(p$data$feature), 3)
    expect_equal(
        vapply(p$layers, function(x) is(x$geom), character(1)),
        paste0("Geom", c("Col", "Label", "Label"))
    )

    p <- plotPie(df, "feature", total_geom = "none")
    expect_equal(
        vapply(p$layers, function(x) is(x$geom), character(1)),
        paste0("Geom", c("Col", "Label"))
    )

    p <- plotPie(df, "feature", cat_geom = "text")
    expect_equal(
        vapply(p$layers, function(x) is(x$geom), character(1)),
        paste0("Geom", c("Col", "Text", "Label"))
    )

})

test_that(".plotDoublePie creates the expected data structures", {

    p <- plotPie(df, "feature", "TF1")
    expect_equal(dim(p$data), c(9, 12))
    expect_equal(
        colnames(p$data),
        c(
            "feature", "TF1", "value", "p", "label_radians", "N", "n", "lab",
            "r", "x", "lab_x", "lab_y"
        )
    )
    expect_equal(
        vapply(p$layers, function(x) is(x$geom), character(1)),
        paste0("Geom", c("ArcBar", "Label", "Label"))
    )
    expect_equal(
        p$labels[c("x", "y", "fill", "r", "label")],
        list(
            x = "TF1", y = "y", fill = "feature", r = "width * r",
            label = "lab"
        )
    )

    p <- plotPie(df, "feature", "TF1", total_geom = "none", cat_geom = "text")
    expect_equal(
        vapply(p$layers, function(x) is(x$geom), character(1)),
        paste0("Geom", c("ArcBar", "Text"))
    )

})

test_that(".plotTriplePie creates the expected data structures", {

    p <- plotPie(df, "feature", "TF1", "TF2")
    expect_equal(dim(p$data), c(27, 14))
    expect_equal(
        colnames(p$data),
        c(
            "feature", "TF1", "TF2", "value", "p", "label_radians", "N", "n",
            "lab", "r", "x", "y", "lab_x", "lab_y"
        )
    )
    expect_equal(
        vapply(p$layers, function(x) is(x$geom), character(1)),
        paste0("Geom", c("ArcBar", "Label", "Label"))
    )
    expect_equal(
        p$labels[c("x", "y", "fill", "r", "label")],
        list(
            x = "TF1", y = "TF2", fill = "feature", r = "width * r", label = "lab"
        )
    )

    p <- plotPie(df, "feature", "TF1", "TF2", total_geom = "none", cat_geom = "text")
    expect_equal(
        vapply(p$layers, function(x) is(x$geom), character(1)),
        paste0("Geom", c("ArcBar", "Text"))
    )
})

test_that("Scaling by columns works as expected", {

    gr <- ex_prom
    mcols(gr) <- df[seq_along(gr),]
    p <- plotPie(gr, fill = "feature")
    expect_equal(sum(p$data$n), length(gr))
    p <- plotPie(gr, fill = "feature", scale_by = "width")
    expect_equal(sum(p$data$n), sum(width(gr) / 1e3))
    df$scale <- 0.5
    p <- plotPie(df, fill = "feature", x = "TF1", scale_by = "scale")
    expect_equal(sum(p$data$value), nrow(df) / 2)
    df$scale <- "a"
    expect_error(plotPie(df, fill = "feature", x = "TF1", scale_by = "scale"))
})

test_that("DataFrame objects work as expected", {
    p <- plotPie(DataFrame(df), fill = "feature")
    expect_true(is(p, "gg"))
    expect_true(is.factor(p$data$feature))
    expect_equal(length(p$data$feature), 3)
    expect_equal(
        vapply(p$layers, function(x) is(x$geom), character(1)),
        paste0("Geom", c("Col", "Label", "Label"))
    )
})
