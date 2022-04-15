bw <- system.file("tests", "test.bw", package = "rtracklayer")
gr <- GRanges("chr2:500")

# Stop testing this to try get build below 10 mins
test_that("Paths behave correctly", {
  expect_s4_class(getProfileData(bw, gr, upstream = 10, bins = 10), "GRanges")
  expect_s4_class(getProfileData(c(bw, bw), gr, upstream = 10, bins = 10), "GRangesList")
})

