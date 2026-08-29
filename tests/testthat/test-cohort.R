test_that("cohort manifests resolve relative paths", {
  cohort <- read_ichor_cohort(fixture_path("example-manifest.csv"), "hg38")

  expect_s3_class(cohort, "ichor_cohort")
  expect_equal(names(cohort$samples), c("example-a", "example-b"))
  expect_equal(cohort$metadata$condition, c("A", "B"))
})

test_that("cohort matrices preserve genomic order and weighted values", {
  cohort <- read_ichor_cohort(fixture_path("example-manifest.csv"), "hg38")
  mat <- ichor_matrix(cohort, bin_size = 1e6, value = "corrected_copy_number")

  expect_s3_class(mat, "ichor_matrix")
  expect_equal(dim(mat$values)[1], 2)
  expect_equal(mat$bins$chr[1:4], rep("1", 4))
  expect_equal(unname(mat$values["example-a", 1:4]), c(2, 2, 3, 3))
  expect_equal(unname(mat$values["example-b", 1:4]), rep(2, 4))
})

test_that("call matrices use documented scores", {
  cohort <- read_ichor_cohort(fixture_path("example-manifest.csv"), "hg38")
  mat <- ichor_matrix(cohort, bin_size = 1e6, value = "call")

  expect_equal(unname(mat$values["example-a", 1:4]), c(0, 0, 1, 1))
})

test_that("heatmap adapter returns a ComplexHeatmap object", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  cohort <- read_ichor_cohort(fixture_path("example-manifest.csv"), "hg38")
  mat <- ichor_matrix(cohort, value = "call")

  expect_s4_class(plot_ichor_heatmap(mat, annotation_columns = "condition"), "Heatmap")
})
