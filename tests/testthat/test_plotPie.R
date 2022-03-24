set.seed(200)
df <- data.frame(
  feature = sample(c("Promoter", "Enhancer", "Intergenic"), 200, replace = TRUE),
  TF1 = sample(c("Up", "Down", "Unchanged"), 200, replace = TRUE),
  TF2 = sample(c("Up", "Down", "Unchanged"), 200, replace = TRUE)
)


test_that("plotPie Errors where expected", {

  expect_error(plotPie(NULL))
  expect_error(plotPie(df), "The initial category must be defined as")
  expect_error(plotPie(df, fill = ""))
  expect_error(plotPie(df, fill = "feature", x = ""))
  expect_error(plotPie(df, fill = "feature", x = "TF1", y = ""))

})

test_that(".plotSinglePie creates expected data structures", {

  p <- plotPie(df, "feature")
  expect_true(is(p, "gg"))
  expect_true(is.factor(p$data$feature))
  expect_equal(length(p$data$feature), 3)
  expect_equal(sum(vapply(p$layers, is, TRUE, "LayerInstance")), 4)

  p <- plotPie(df, "feature", show_total = FALSE)
  expect_equal(sum(vapply(p$layers, is, TRUE, "LayerInstance")), 3)

  p <- plotPie(df, "feature", show_category = FALSE)
  expect_equal(sum(vapply(p$layers, is, TRUE, "LayerInstance")), 2)

})

test_that(".plotDoublePie creates the expected data structures", {

  p <- plotPie(df, "feature", "TF1")
  expect_equal(dim(p$data), c(3, 7))
  expect_equal(
    colnames(p$data),
    c("TF1", "N", "r", "Enhancer", "Intergenic", "Promoter", "x")
  )
  expect_equal(sum(vapply(p$layers, is, TRUE, "LayerInstance")), 2)
  expect_equal(
    p$labels[c("x", "y", "fill", "r", "label")],
    list(
      x = "TF1", y = "y", fill = "feature", r = "width * r",
      label = "comma(N, 1)"
    )
  )

  p <- plotPie(df, "feature", "TF1", show_total = FALSE)
  expect_equal(sum(vapply(p$layers, is, TRUE, "LayerInstance")), 1)

})

test_that(".plotDoublePie creates the expected data structures", {

  p <- plotPie(df, "feature", "TF1", "TF2")
  expect_equal(sum(vapply(p$layers, is, TRUE, "LayerInstance")), 2)
  expect_equal(dim(p$data), c(9, 9))
  expect_equal(
    colnames(p$data),
    c("TF1", "TF2", "N", "r", "Enhancer", "Intergenic", "Promoter", "x", "y")
  )
  expect_equal(
    p$labels[c("x", "y", "fill", "r", "label")],
    list(
      x = "TF1", y = "TF2", fill = "feature", r = "width * r",
      label = "comma(N, 1)"
    )
  )

  p <- plotPie(df, "feature", "TF1", "TF2", show_total = FALSE)
  expect_equal(sum(vapply(p$layers, is, TRUE, "LayerInstance")), 1)

})
