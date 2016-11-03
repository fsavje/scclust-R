
library(Rscclust)
context("Appveyor")

test_that("Appveyor skips as expected", {
  skip_on_travis()
  skip_on_appveyor()
  expect_silent(stop("hej"))
})
