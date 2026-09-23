test_that("sample objects validate and print", {
  x <- example_sample("a")

  expect_s3_class(x, "ichor_sample")
  expect_equal(x$sample_id, "example-a")
  expect_equal(x$genome_build, "hg38")
  expect_equal(nrow(x$segments), 8)
  expect_invisible(validate_ichor_sample(x))
  expect_output(print(x), "12")
})

test_that("tumor fractions use fractional units", {
  x <- example_sample("a")
  x$params$tumor_fraction <- 12.5
  expect_error(validate_ichor_sample(x), "fraction", class = "seena_validation_error")
})

test_that("input provenance is available", {
  p <- ichor_provenance(example_sample("a"))
  expect_equal(sort(p$type), c("cna_seg", "params", "seg"))
  expect_true(all(nchar(p$md5) == 32))
})
