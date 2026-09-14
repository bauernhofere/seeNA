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

test_that("all profile plotting paths validate point size consistently", {
  a <- example_sample("a")
  for (bad in list(-1, NA_real_, Inf, c(0, 1), "small")) {
    expect_error(plot_ichor_profile(a, point_size = bad), "point_size")
    expect_error(plot_ichor_compare(a, point_size = bad), "point_size")
    expect_error(plot_ichor_region(a, "chr1", point_size = bad), "point_size")
  }
  expect_s3_class(plot_ichor_profile(a, point_size = 0), "ggplot")
})

test_that("sample identifiers are shown by default and can be hidden", {
  a <- example_sample("a")
  b <- example_sample("b")

  shown <- plot_ichor_profile(a)
  expect_identical(shown$labels$title, "example-a")
  hidden <- plot_ichor_profile(a, show_sample_id = FALSE)
  expect_null(hidden$labels$title)
  expect_false(grepl("example-a", hidden$labels$subtitle, fixed = TRUE))

  compare <- plot_ichor_compare(list(a, b))
  expect_match(compare$labels$subtitle, "example-a")
  built <- ggplot2::ggplot_build(compare)
  expect_identical(ggplot2::get_guide_data(built, "colour")$.label, c("example-a", "example-b"))

  anonymous <- plot_ichor_compare(list(a, b), show_sample_id = FALSE)
  expect_identical(anonymous$labels$subtitle, "hg38")
  expect_null(ggplot2::get_guide_data(ggplot2::ggplot_build(anonymous), "colour"))
  region <- plot_ichor_region(list(a, b), "chr1:1-4000000", show_sample_id = FALSE)
  expect_false(grepl("example", region$labels$subtitle, fixed = TRUE))
  expect_error(plot_ichor_profile(a, show_sample_id = NA), "TRUE or FALSE")
})
