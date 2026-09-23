segment_cohort <- function() {
  read_ichor_cohort(fixture_path("example-manifest.csv"), "hg38")
}

test_that("segment matrices use exported medians and leave bins and provenance untouched", {
  c <- segment_cohort()
  original <- c
  m <- ichor_matrix(c, value = "segment_median", chromosomes = "1")
  expect_equal(unname(m$values[1, 1:4]), c(-0.005, -0.005, 0.58, 0.58))
  expect_equal(unname(m$values[2, 1:4]), rep(0.03, 4))
  expect_equal(unname(m$coverage[, 1:4]), matrix(1, 2, 4))
  expect_false(any(m$mixed))
  expect_identical(m$value, "segment_median")
  expect_identical(m$aggregation, "base_pair_mean")
  expect_null(m$call_column)
  expect_identical(m$provenance, ichor_provenance(c))
  expect_identical(c, original)
  expect_false(identical(m$values, ichor_matrix(c, value = "logR", chromosomes = "1")$values))
  expect_invisible(validate_ichor_matrix(m))
})

test_that("segment boundaries use inclusive overlap-weighted means", {
  c <- segment_cohort()
  s <- c$samples[[1]]$segments[1:2, ]
  s$start <- c(1, 1250001); s$end <- c(1250000, 3000000)
  s$median <- c(-0.4, 0.8)
  c$samples[[1]]$segments <- s
  m <- ichor_matrix(c, value = "segment_median", chromosomes = "1")
  expect_equal(unname(m$values[1, 1:3]), c(-0.4, 0.5, 0.8))
  expect_true(is.na(m$values[1, 4]))
  expect_equal(m$coverage[1, 4], 0, ignore_attr = TRUE)
  fine <- ichor_matrix(c, value = "segment_median", chromosomes = "1", bin_size = 5e5)
  expect_equal(unname(fine$values[1, 1:6]), c(-0.4, -0.4, 0.2, 0.8, 0.8, 0.8))
})

test_that("missing medians and segment gaps obey the same coverage threshold", {
  c <- segment_cohort()
  s <- c$samples[[1]]$segments[1:2, ]
  s$start <- c(1, 1250001); s$end <- c(1250000, 3000000)
  s$median <- c(NA_real_, 0.8)
  c$samples[[1]]$segments <- s
  full <- ichor_matrix(c, value = "segment_median", chromosomes = "1")
  expect_equal(unname(full$coverage[1, 1:4]), c(0, 0.75, 1, 0))
  expect_true(all(is.na(full$values[1, c(1, 2, 4)])))
  partial <- ichor_matrix(c, value = "segment_median", chromosomes = "1", min_coverage = 0.75)
  expect_equal(partial$values[1, 2], 0.8, ignore_attr = TRUE)
  zero <- ichor_matrix(c, value = "segment_median", chromosomes = "1", min_coverage = 0)
  expect_true(all(is.na(zero$values[1, c(1, 4)])))
  c$samples[[1]]$segments$median[] <- NA_real_
  empty <- ichor_matrix(c, value = "segment_median", chromosomes = "1")
  expect_true(all(is.na(empty$values[1, ])))
  expect_true(all(empty$coverage[1, ] == 0))
})

test_that("internal segment gaps are unsupported and zero medians remain observed", {
  c <- segment_cohort()
  s <- c$samples[[1]]$segments[1:2, ]
  s$start <- c(1, 1750001); s$end <- c(1250000, 3000000)
  s$median <- c(0, 0.8)
  c$samples[[1]]$segments <- s
  c$samples[[1]]$params <- NULL
  m <- ichor_matrix(c, value = "segment_median", chromosomes = "1")
  expect_equal(unname(m$values[1, 1]), 0)
  expect_true(is.na(m$values[1, 2]))
  expect_equal(unname(m$coverage[1, 2]), 0.5)
  partial <- ichor_matrix(c, value = "segment_median", chromosomes = "1", min_coverage = 0.5)
  expect_equal(unname(partial$values[1, 2]), 0.4)
})

test_that("segment support is not an undocumented bin-logR mask", {
  c <- segment_cohort()
  c$samples[[1]]$bins <- c$samples[[1]]$bins[-2, ]
  c$samples[[1]]$bins$logR[1] <- NA_real_
  segment <- ichor_matrix(c, value = "segment_median", chromosomes = "1")
  bins <- ichor_matrix(c, value = "logR", chromosomes = "1")
  expect_equal(unname(segment$values[1, 1:2]), c(-0.005, -0.005))
  expect_equal(unname(segment$coverage[1, 1:2]), c(1, 1))
  expect_true(all(is.na(bins$values[1, 1:2])))
})

test_that("missing segment files error, but absent chromosomes remain missing", {
  c <- segment_cohort()
  c$samples[[2]]$segments <- NULL
  expect_error(ichor_matrix(c, value = "segment_median"), "example-b.*segment")
  expect_s3_class(ichor_matrix(c, value = "logR"), "ichor_matrix")
  c <- segment_cohort()
  c$samples[[1]]$segments <- c$samples[[1]]$segments[c$samples[[1]]$segments$chr == "1", ]
  m <- ichor_matrix(c, value = "segment_median", chromosomes = "X")
  expect_true(all(is.na(m$values[1, ])))
  expect_true(all(m$coverage[1, ] == 0))
})

test_that("short terminal target bins use actual chromosome widths", {
  c <- segment_cohort()
  limit <- ichor_genome_layout("hg38", "19")$length
  for (i in seq_along(c$samples)) {
    s <- c$samples[[i]]$segments[1, ]
    s$chr <- "19"; s$start <- 58000001; s$end <- limit; s$median <- -0.25
    c$samples[[i]]$segments <- s
  }
  m <- ichor_matrix(c, value = "segment_median", chromosomes = "19")
  expect_equal(unname(m$values[, ncol(m$values)]), c(-0.25, -0.25))
  expect_equal(unname(m$coverage[, ncol(m$values)]), c(1, 1))
  expect_true(all(is.na(m$values[, -ncol(m$values)])))
})

test_that("segment matrices reject invalid sources and respect allocation limits", {
  c <- segment_cohort()
  expect_error(ichor_matrix(c, value = "segment_median", max_cells = 1), "max_cells")
  bad <- c; bad$samples[[1]]$segments$median[1] <- Inf
  expect_error(ichor_matrix(bad, value = "segment_median"), "invalid median")
  bad <- c; bad$samples[[1]]$segments$start[2] <- 2000000
  expect_error(ichor_matrix(bad, value = "segment_median"), "overlapping")
  m <- ichor_matrix(c, value = "segment_median", chromosomes = "1")
  m$mixed[1, 1] <- TRUE
  expect_error(validate_ichor_matrix(m), "Mixed-call flags")
})

test_that("segment medians render through the ComplexHeatmap adapter", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  c <- segment_cohort()
  m <- ichor_matrix(c, value = "segment_median", chromosomes = "1")
  h <- plot_ichor_heatmap(m, annotation_columns = "condition")
  expect_identical(h@matrix, m$values)
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  tryCatch(expect_s4_class(ComplexHeatmap::draw(h), "HeatmapList"), finally = grDevices::dev.off())
  expect_gt(file.info(f)$size, 0)
})
