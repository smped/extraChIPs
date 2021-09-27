dat <- cbind(
  matrix(rnorm(1000), nrow=100, ncol=5),
  matrix(rnorm(1000, .1, .7), nrow=100, ncol=5)
)
colnames(dat) <- c(paste0("A", 1:5), paste0("B", 1:5))
df <- DataFrame(group = rep(c("A", "B"), each = 5))
se <- SummarizedExperiment(
  assays = SimpleList(counts = dat),
  colData = df
)

test_that("factor is handled correctly", {

  expect_error(addQSmooth(se, "counts", "missing"))
  expect_error(addQSmooth(se, "counts"))
  expect_s4_class(addQSmooth(se, factor = "group"), "SummarizedExperiment")

})
