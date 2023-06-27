set.seed(100)
gr1 <- GRanges(paste0("chr1:", seq(10, 150, by = 10)))
width(gr1) <- 5
gr1$logFC <- rnorm(length(gr1))
gr1$status <- dplyr::case_when(
  gr1$logFC > 0.5 ~ "Increased",
  gr1$logFC < -0.5 ~ "Decreased",
  TRUE ~ "Unchanged"
)
gr2 <-  GRanges(paste0("chr1:", seq(51, 250, by = 15)))
width(gr2) <- 4
gr2$logFC <- rnorm(length(gr2))
gr2$status <- dplyr::case_when(
  gr2$logFC > 0.5 ~ "Increased",
  gr2$logFC < -0.5 ~ "Decreased",
  TRUE ~ "Unchanged"
)
grl <- GRangesList(TF1 = gr1, TF2 = gr2)

test_that(
  "plotPairwise works correctly", {
    p <- plotPairwise(grl, var = "logFC", colour = "status")
    expect_true(is(p, "ggside"))

    p <- plotPairwise(grl, var = "logFC")
    expect_true(is(p, "ggside"))

    p <- plotPairwise(grl, var = "logFC", xside = "none", yside = "none")
    expect_true(is(p, "gg"))

    gr1$status <- as.factor(gr1$status)
    gr2$status <- as.factor(gr2$status)
    grl <- GRangesList(TF1 = gr1, TF2 = gr2)
    p <- plotPairwise(grl, var = "logFC", colour = "status", xside = "density", yside = "density")
    expect_true(is(p, "ggside"))
    p <- plotPairwise(grl, var = "logFC", colour = "status", xside = "violin", yside = "violin")
    expect_true(is(p, "ggside"))
  }
)

test_that(
  "plotPairwise errors as expected", {
    expect_error(
      plotPairwise(grl, var = "logFC", colour = "logFC")
    )
    expect_error(
      plotPairwise(grl, var = "status")
    )
  }
)
