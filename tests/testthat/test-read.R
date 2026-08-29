test_that("prefixed cna columns are normalized", {
  d <- read_ichor_cna(fixture_path("example-a.cna.seg"))

  expect_equal(attr(d, "sample_id"), "example-a")
  expect_equal(nrow(d), 12)
  expect_true(all(c("chr", "start", "end", "logR", "copy_number",
                    "corrected_copy_number", "corrected_call") %in% names(d)))
  expect_equal(unique(d$chr), c("1", "2", "X"))
  expect_equal(d$corrected_copy_number[3], 3)
})

test_that("chr prefixes and companion files are normalized", {
  cna <- read_ichor_cna(fixture_path("example-b.cna.seg"))
  seg <- read_ichor_segments(fixture_path("example-b.seg"))
  params <- read_ichor_params(fixture_path("example-b.params.txt"))

  expect_equal(unique(cna$chr), c("1", "2", "X"))
  expect_equal(unique(seg$chr), c("1", "2", "X"))
  expect_equal(params$sample_id, "example-b")
  expect_equal(params$tumor_fraction, 0.072)
  expect_equal(params$ploidy, 1.98)
})

test_that("missing required columns fail early", {
  path <- tempfile(fileext = ".cna.seg")
  writeLines("chr\tstart\tend\n1\t1\t100\n", path)
  expect_error(read_ichor_cna(path), class = "ichorviz_schema_error")
})
