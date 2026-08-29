test_that("profile and region functions return ggplots", {
  a <- example_sample("a")
  b <- example_sample("b")

  expect_s3_class(plot_ichor_profile(a), "ggplot")
  expect_s3_class(plot_ichor_compare(list(a, b)), "ggplot")
  expect_s3_class(plot_ichor_region(list(a, b), "chr1:1-4000000"), "ggplot")
})

test_that("comparison validates sample colors and regions", {
  a <- example_sample("a")
  b <- example_sample("b")

  expect_error(plot_ichor_compare(list(a, b), colors = c("example-a" = "red")),
               "named for every sample")
  expect_error(plot_ichor_region(a, "chr8:1-4000000"), "No bins overlap")
})
