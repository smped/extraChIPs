set.seed(101)
df <- data.frame(logFC = rnorm(20), p = rbeta(20, shape1 = 1, shape2 = 20))
df$FDR <- p.adjust(df$p, "fdr")
df$FDR[1] <- NA

gr1 <- GRanges(paste0("chr1:", seq(10, 150, by = 10)))
gr1$logFC <- rnorm(length(gr1))
gr1$FDR <- rbeta(length(gr1), 1, 8)
gr2 <-  GRanges(paste0("chr1:", seq(51, 250, by = 15)))
gr2$logFC <- rnorm(length(gr2))
gr2$FDR <- rbeta(length(gr2), 1, 8)
grl <- GRangesList(TF1 = gr1, TF2 = gr2)

test_that("addDiffStatus performs correctly", {
  nm <- names(df)
  df <- addDiffStatus(df)
  lv <- c("Unchanged", "Decreased", "Increased", "Undetected")
  expect_equal(colnames(df), c(nm, "status"))
  expect_equal(levels(df$status), lv)
  expect_equal(
    levels( addDiffStatus(df, missing = NA_character_)$status), lv[1:3]
  )
  expect_true(is(df, "data.frame"))
  grl <- addDiffStatus(grl)
  expect_true(is(grl, "GRangesList"))
  expect_equal(colnames(mcols(grl[[1]])),  c("logFC", "FDR", "status"))
  expect_equal(
    levels(addDiffStatus(df[-1,], drop = TRUE)$status), lv[1:3]
  )
})
